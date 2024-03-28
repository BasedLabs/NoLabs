from typing import List

from pydantic import BaseModel
from rdkit import Chem
import openeye.oechem as oechem
import openeye.oeomega as oeomega

from copy import deepcopy

from typing_extensions import Literal

from dockstream.utils.dockstream_exceptions import LigandPreparationFailed

from dockstream.core.ligand_preparator import LigandPreparator, _LE
from dockstream.core.RDkit.RDkit_ligand_preparator import RDkitLigandPreparator
from dockstream.utils.translations.molecule_translator import MoleculeTranslator
from dockstream.utils.translations.translation import RDkitMolToOpenEyeMol
from dockstream.utils.enums.OpenEye_enums import OpenEyeLigandPreparationEnum
from dockstream.core.ligand.ligand import Ligand

_LP = OpenEyeLigandPreparationEnum()


class OpenEyeLigandPreparator(LigandPreparator, BaseModel):

    type: Literal["OpenEye"] = "OpenEye"

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)

    def _initialize_ligands(self):
        super()._initialize_ligands()

    def _load_references(self):
        references = []

        for path in self.align.reference_paths:
            mol_supplier = oechem.oemolistream()

            # set the provided format
            ref_format = self.align.reference_format.upper()
            if ref_format == _LP.ALIGN_REFERENCE_FORMAT_SDF:
                mol_supplier.SetFormat(oechem.OEFormat_SDF)
            elif ref_format == _LP.ALIGN_REFERENCE_FORMAT_PDB:
                mol_supplier.SetFormat(oechem.OEFormat_PDB)
            else:
                raise LigandPreparationFailed("Specified format not supported!")

            if mol_supplier.open(path):
                for mol in mol_supplier.GetOEMols():
                    references.append(oechem.OEMol(mol))
            else:
                oechem.OEThrow.Fatal("Unable to create specified output file.")
        if len(references) == 0:
            raise LigandPreparationFailed("No reference molecules could be loaded at path(s) specified.")
        self._references = references
        self._logger.log(f"Stored {len(references)} reference molecules.", _LE.DEBUG)

    def _get_RDkit_aligner(self, conf, ligands):
        return RDkitLigandPreparator(ligands=ligands, **conf)

    def _smiles_to_molecules(self, ligands: List[Ligand]) -> List[Ligand]:
        for lig in ligands:
            lig_molecule = oechem.OEMol()
            oechem.OESmilesToMol(lig_molecule, lig.get_smile())
            lig.set_molecule(lig_molecule)
            lig.set_mol_type(_LP.TYPE_OPENEYE)
        return ligands

    def generate3Dcoordinates(self):
        """Method to generate 3D coordinates, in case the molecules have been built from SMILES."""

        for lig in self.ligands:
            lig.set_molecule(None)
            lig.set_mol_type(None)
        ligand_list = self._smiles_to_molecules(deepcopy(self.ligands))

        failed = 0
        succeeded = 0
        builder = oeomega.OEConformerBuilder()
        for idx, ligand in enumerate(ligand_list):
            inp_mol = ligand.get_molecule()
            if inp_mol is None:
                continue
            return_code = builder.Build(inp_mol)
            if return_code != oeomega.OEOmegaReturnCode_Success:
                failed += 1
                self._logger.log(f"The 3D coordinate generation of molecule {ligand.get_ligand_number()} (smile: {ligand.get_smile()}) failed (oeomega return code={return_code}).",
                                 _LE.DEBUG)
                continue
            self.ligands[idx] = Ligand(smile=ligand.get_smile(),
                                       original_smile=ligand.get_original_smile(),
                                       ligand_number=ligand.get_ligand_number(),
                                       enumeration=ligand.get_enumeration(),
                                       molecule=oechem.OEMol(inp_mol),
                                       mol_type=_LP.TYPE_OPENEYE,
                                       name=ligand.get_name())
            succeeded += 1

        if failed > 0:
            self._logger.log(f"Of {len(self.ligands)}, {failed} could not be embedded.", _LE.WARNING)
        self._logger.log(f"In total, {succeeded} ligands were successfully embedded (oeomega).", _LE.DEBUG)

    def align_ligands(self):
        if self.align.mode != _LP.ALIGN_MODE_INTERNAL:
            raise LigandPreparationFailed("Only internal alignment supported at the moment.")
        if self._references is None:
            raise LigandPreparationFailed("No reference molecule has been found.")

        # use the general, internal alignment technique
        # ---------
        # 1) translate the ligands from openeye to rdkit and do not use "bySMILES" method, as
        #    coordinates would be lost
        mol_trans = MoleculeTranslator(self.ligands)
        ligands_rdkit = mol_trans.get_as_rdkit()
        self._logger.log(f"Align: Of {len(self.ligands)}, {len(ligands_rdkit)} were translated to RDkit molecules.",
                         _LE.DEBUG)

        # 2) do the alignment to a reference molecule; also disable RDkit logger
        ligands_rdkit = self._align_ligands_with_RDkit_preparator(ligands_rdkit)

        # 3) translate ligands back and update internal collection
        mol_trans = MoleculeTranslator(ligands_rdkit)
        translated_mols = mol_trans.get_as_openeye()
        for lig, translated_mol in zip(self.ligands, translated_mols):
            lig.set_molecule(translated_mol.get_molecule())

    def write_ligands(self, path, format):
        ofs = oechem.oemolostream()
        format = format.upper()

        ligands_copy = [deepcopy(lig) for lig in self.ligands]

        # check and specify format of file
        if format == _LP.OUTPUT_FORMAT_SDF:
            ofs.SetFormat(oechem.OEFormat_SDF)
        elif format == _LP.OUTPUT_FORMAT_MOL2:
            ofs.SetFormat(oechem.OEFormat_MOL2)
        else:
            raise LigandPreparationFailed("Specified output format unknown.")

        if ofs.open(path):
            for lig in ligands_copy:
                lig.add_tags_to_molecule()
                if lig.get_molecule() is not None:
                    mol = deepcopy(lig.get_molecule())
                    mol.SetTitle(lig.get_identifier())
                    oechem.OEWriteMolecule(ofs, mol)
        else:
            oechem.OEThrow.Fatal("Unable to create specified output file.")
        ofs.close()
        self._logger.log(f"Wrote {len(self.ligands)} molecules to file {path} (format: {format}).", _LE.DEBUG)

    def _make_ligands_from_molecules(self, ligands):
        buffer = []
        if isinstance(ligands[0], Chem.Mol):
            ligands = [RDkitMolToOpenEyeMol(mol, bySMILES=False) for mol in ligands]
        for index_mol, mol in enumerate(ligands):
            buffer.append(Ligand(smile=oechem.OEMolToSmiles(mol),
                                 ligand_number=index_mol,
                                 enumeration=0,
                                 molecule=mol,
                                 mol_type=_LP.TYPE_OPENEYE))
        self.ligands = buffer
