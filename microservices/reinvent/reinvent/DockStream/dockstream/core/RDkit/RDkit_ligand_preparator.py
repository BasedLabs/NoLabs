import os
from copy import deepcopy

from typing import Optional, List
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from typing_extensions import Literal

from dockstream.utils.dockstream_exceptions import LigandPreparationFailed

from dockstream.core.ligand_preparator import LigandPreparator, _LE
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum
from dockstream.utils.smiles import to_mol
from dockstream.core.ligand.ligand import Ligand

_LP = RDkitLigandPreparationEnum()

# This class is inspired / based on code written by
# Peter Schmidtke (https://github.com/Discngine/rdkit_tethered_minimization), which is accessible under the MIT license
# and the blogspot entry called "more-on-constrained-embedding" from the RDkit guys.


class ParametersCoordinateGeneration(BaseModel):
    method: str = "UFF"
    maximum_iterations: Optional[int] = 600


class RDkitLigandPreparatorParameters(BaseModel):
    protonate: Optional[bool] = True
    coordinate_generation: ParametersCoordinateGeneration = ParametersCoordinateGeneration()


class RDkitLigandPreparator(LigandPreparator, BaseModel):
    """Class that deals with all the preparatory steps needed before actual docking using "rDock" can commence."""

    type: Literal["RDkit"] = "RDkit"
    parameters: RDkitLigandPreparatorParameters = RDkitLigandPreparatorParameters()

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)

    def _load_references(self):
        references = []
        ref_format = self.align.reference_format.upper()
        for path in self.align.reference_paths:
            if ref_format == _LP.ALIGN_REFERENCE_FORMAT_PDB:
                ref_mol = Chem.MolFromPDBFile(path, sanitize=True)
                ref_mol.SetProp("_Name", os.path.basename(path))
                references.append(ref_mol)
            elif ref_format == _LP.ALIGN_REFERENCE_FORMAT_SDF:
                mol_supplier = Chem.SDMolSupplier(path)
                for mol in mol_supplier:
                    references.append(mol)
            else:
                raise IOError("Specified format not supported.")
        if len(references) == 0:
            raise LigandPreparationFailed("No reference molecules could be loaded with path(s) specified.")
        self._references = references
        self._logger.log(f"Stored {len(references)} reference molecules.", _LE.DEBUG)

    def _smiles_to_molecules(self, ligands: List[Ligand]) -> List[Ligand]:
        for lig in ligands:
            mol = to_mol(lig.get_smile())
            lig.set_molecule(mol)
            lig.set_mol_type(_LP.TYPE_RDKIT)
        return ligands

    def generate3Dcoordinates(self, converged_only=False):
        """Method to generate 3D coordinates, in case the molecules have been built from SMILES."""

        for lig in self.ligands:
            lig.set_molecule(None)
            lig.set_mol_type(None)
        ligand_list = self._smiles_to_molecules(deepcopy(self.ligands))

        failed = 0
        succeeded = 0
        for idx, lig_obj in enumerate(ligand_list):
            ligand = lig_obj.get_molecule()
            if ligand is None:
                continue

            # note, that parameter "useRandomCoords" needs to be "True", which is often required for larger molecules
            # as the embedding sometimes fails
            embed_code = AllChem.EmbedMolecule(ligand, randomSeed=42, useRandomCoords=True)

            # while MMFF sometimes gives better geometries, UFF has a wider range of parameters and thus will fail less
            # often and is also much quicker
            status = 0
            if self.parameters.coordinate_generation.method == _LP.EP_PARAMS_COORDGEN_UFF:
                # check, if embedding worked
                if embed_code != -1:
                    status = AllChem.UFFOptimizeMolecule(ligand, maxIters=self.parameters.coordinate_generation.maximum_iterations)
                    if status == 1:
                        self._logger.log(f"The 3D coordinate generation of molecule number {lig_obj.get_ligand_number()} (smile: {lig_obj.get_smile()}) did not converge in time - try increasing the number of maximum iterations.",
                                         _LE.DEBUG)
                        failed += 1
                        if converged_only:
                            continue
                else:
                    self._logger.log(f"Could not embed molecule number {lig_obj.get_ligand_number()} (smile: {lig_obj.get_smile()}) - no 3D coordinates generated.",
                                     _LE.DEBUG)
                    failed += 1
                    continue
            else:
                raise LigandPreparationFailed("Coordination generation method %s is not supported." % self.parameters.coordinate_generation.method)

            # add hydrogens to the molecule
            if self.parameters.protonate:
                ligand = Chem.AddHs(ligand, addCoords=True)

            self.ligands[idx] = Ligand(smile=lig_obj.get_smile(),
                                       original_smile=lig_obj.get_original_smile(),
                                       ligand_number=lig_obj.get_ligand_number(),
                                       enumeration=lig_obj.get_enumeration(),
                                       molecule=ligand,
                                       mol_type=lig_obj.get_mol_type(),
                                       name=lig_obj.get_name())
            succeeded += 1

        if failed > 0:
            self._logger.log(f"Of {len(self.ligands)}, {failed} could not be embedded.",
                             _LE.WARNING)
        self._logger.log(f"In total, {succeeded} ligands were successfully embedded (RDkit).", _LE.DEBUG)

    def align_ligands(self):
        """This method loops over the molecules stored and structurally aligns them to a reference molecule. If
           a list longer is provided, the reference molecule with the largest structural overlap is chosen."""

        molecules_container = []
        failed = 0

        # make copies of the ligands in case something goes wrong; also remove hydrogens in case they are present
        # as they hinder core calculation and are added afterwards
        #ligands = [Chem.RemoveHs(Chem.Mol(deepcopy(molecule))) for molecule in [lig.get_molecule() for lig in self.ligands]]
        references = [Chem.RemoveHs(Chem.Mol(reference)) for reference in deepcopy(self._references)]
        for lig_obj in self.ligands:
            success = False
            if lig_obj.get_molecule() is None:
                continue
            ligand = Chem.RemoveHs(Chem.Mol(deepcopy(lig_obj.get_molecule())))

            # find the maximum common substructure for the current molecule and the reference molecule with the largest
            # topological overlap
            most_mcs_atoms = 0
            for cur_reference in references:
                # in this implementation, a threshold of > 0.5 will always lead to a MCS only to be returned if both
                # molecules have a FMC, but for future implementations this might change
                cur_mcs = rdFMCS.FindMCS([cur_reference, ligand],
                                         threshold=1.0,
                                         completeRingsOnly=self.align.complete_rings_only)

                # check, if this reference has a larger overlap (in atoms for the MCS) than the best hit so far
                if cur_mcs.numAtoms > most_mcs_atoms:
                    reference = cur_reference
                    mcs = cur_mcs
                    most_mcs_atoms = cur_mcs.numAtoms

            # check, if a substructure could be found and proceed accordingly; note, that an empty SMARTS string is
            # entirely meaningless
            if most_mcs_atoms > 0 and mcs.smartsString and len(mcs.smartsString) > 0:
                # generate a molecule from the MCS SMART to extract that part from the reference later
                mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

                if mcs_mol:
                    # generate a copy of the molecule as the coordinates will be updated now
                    new_molecule = Chem.Mol(ligand)

                    # get the matches for both the reference molecule and the current ligand; this is necessary, as
                    # in some cases the MCS will not be recognised otherwise as atom specification might differ
                    reference_match = reference.GetSubstructMatch(mcs_mol)
                    ligand_match = new_molecule.GetSubstructMatch(mcs_mol)

                    # align the current ligand to the reference structure
                    Chem.AllChem.AlignMol(new_molecule, reference, atomMap=list(zip(ligand_match,
                                                                                    reference_match)))

                    # generate the ratio of the number of core atoms to the number of all reference atoms
                    match_ratio = float(mcs_mol.GetNumAtoms()) / float(reference.GetNumAtoms())

                    # generate the IDs of the atoms for tethering later
                    tethered_atom_IDs = new_molecule.GetSubstructMatches(mcs_mol)

                    if tethered_atom_IDs and match_ratio >= self.align.minimum_substructure_ratio:
                        # get only the list of atoms from tupel
                        tethered_atom_IDs = tethered_atom_IDs[0]

                        # change indexing from 0-indexed to 1-indexed and generate string for property assignment
                        tethered_atom_IDs = map(lambda x: x + 1, list(tethered_atom_IDs))
                        tethered_atom_IDs = ','.join(str(cur_id) for cur_id in tethered_atom_IDs)

                        # add hydrogens and their coordinates
                        new_molecule = Chem.AddHs(new_molecule,
                                                  addCoords=True)

                        # meta-info: (1) atoms, which were used in the alignment,
                        #            (2) name of the reference used (if set, otherwise it will be an empty string),
                        #            (3) the ratio of the matched atoms on the overall atoms of the reference
                        new_molecule.SetProp(_LP.TAG_ALIGNED_ATOMS, tethered_atom_IDs)
                        new_molecule.SetProp(_LP.TAG_ALIGNED_REFERENCE, reference.GetProp("_Name"))
                        new_molecule.SetProp(_LP.TAG_ALIGNED_ATOMS_RATIO,
                                             str(round(match_ratio, ndigits=3)))
                        molecules_container.append(Ligand(smile=lig_obj.get_smile(),
                                                          original_smile=lig_obj.get_original_smile(),
                                                          ligand_number=lig_obj.get_ligand_number(),
                                                          enumeration=lig_obj.get_enumeration(),
                                                          molecule=new_molecule,
                                                          mol_type=lig_obj.get_mol_type(),
                                                          name=lig_obj.get_name()))
                        success = True

            # handle cases, where the molecule could not be aligned to any reference under the
            # threshold constraints given
            if not success:
                failed += 1
                if self.align.fail_action == _LP.ALIGN_FAIL_KEEP:
                    self._logger.log(f"Ligand {repr(lig_obj)} could not be aligned - keeping old coordinates.",
                                     _LE.DEBUG)
                elif self.align.fail_action == _LP.ALIGN_FAIL_DISCARD:
                    self._logger.log(f"Ligand {repr(lig_obj)} could not be aligned - removing coordinates.",
                                     _LE.DEBUG)
                    lig_obj.set_molecule(None)
                    lig_obj.set_mol_type(None)
                else:
                    raise LigandPreparationFailed("Specified align fail action is not implemented.")
                molecules_container.append(lig_obj)

        # replace the molecules stored internally with the aligned ones
        self.ligands = molecules_container
        self._logger.log(f"Attempted to align {len(molecules_container)} enumerations to a reference molecule - {failed} failed (align fail action={_LP.ALIGN_FAIL_KEEP}).", _LE.DEBUG)

        if self.align.tethering:
            self._tether_ligands()

    def _tether_ligands(self):
        # check, if alignment has been done and do it if not (but give a warning to the user)
        # note, that if properties in RDkit molecules are called that are not there, it will result in a KeyError
        for lig in self.ligands:
            try:
                _ = lig.get_molecule().GetProp(_LP.TAG_ALIGNED_ATOMS)
            except KeyError:
                self._logger.log("Call of _tether_ligands() before molecules are aligned - will try to align them now.", _LE.WARNING)
                super().align_ligands()
                break

        # for rDock, specify the atoms to which the molecule is to be tethered during docking under a
        # hard-coded tag (from source in rDock); here, use the atoms used for alignment earlier
        for lig in self.ligands:
            lig.get_molecule().SetProp(_LP.TAG_RDOCK_TETHERED_ATOMS,
                                       lig.get_molecule().GetProp(_LP.TAG_ALIGNED_ATOMS))

        self._logger.log(f"Set tethering tag (rDock) for {len(self.ligands)} molecules.", _LE.DEBUG)
