from typing import List
from dockstream.utils.translations.translation import RDkitMolToOpenEyeMol, OpenEyeMolToRDkitMol
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.core.ligand.ligand import Ligand


class MoleculeTranslator:

    def __init__(self, molecules: List[Ligand], force_mol_type=False):
        self._LP = LigandPreparationEnum()
        self._known_types = [self._LP.TYPE_RDKIT, self._LP.TYPE_OPENEYE, self._LP.TYPE_CORINA,
                             self._LP.TYPE_GOLD, self._LP.TYPE_LIGPREP, self._LP.TYPE_OMEGA]

        if not isinstance(molecules, list):
            molecules = [molecules]
        self._molecules = molecules
        if force_mol_type is False:
            mol_type = molecules[0].get_mol_type()
        else:
            mol_type = force_mol_type
        if mol_type is None:
            raise Exception("Type None cannot be used in molecule translations.")
        if mol_type not in self._known_types:
            raise ValueError(f"Type {mol_type} not in list of supported types.")
        self._mol_type = mol_type

    def get_as_rdkit(self):
        return self._translate_molecules(self._molecules, from_type=self._mol_type, to_type=self._LP.TYPE_RDKIT)

    def get_as_openeye(self):
        return self._translate_molecules(self._molecules, from_type=self._mol_type, to_type=self._LP.TYPE_OPENEYE)

    def add_molecules(self, molecules: list):
        molecules = [mol.get_clone() for mol in molecules]
        for molecule in molecules:
            mol_type = molecule.get_mol_type()
            self._molecules = self._molecules + self._translate_molecules([molecule], mol_type, self._mol_type)

    def _translate_molecules(self, molecules, from_type, to_type, bySMILES=False) -> list:
        # TODO: cover case, where conformers have been added
        if from_type == to_type or from_type is None:
            return molecules
        else:
            buffer = []
            for mol in molecules:
                if (from_type == self._LP.TYPE_RDKIT or from_type == self._LP.TYPE_CORINA or from_type == self._LP.TYPE_LIGPREP) and \
                to_type == self._LP.TYPE_OPENEYE:
                    buffer.append(Ligand(smile=mol.get_smile(), ligand_number=mol.get_ligand_number(),
                                         original_smile=mol.get_original_smile(),
                                         enumeration=mol.get_enumeration(),
                                         molecule=RDkitMolToOpenEyeMol(mol.get_molecule(), bySMILES=bySMILES),
                                         mol_type=self._LP.TYPE_OPENEYE,
                                         name=mol.get_name()))
                elif from_type == self._LP.TYPE_OPENEYE and \
                (to_type == self._LP.TYPE_RDKIT or to_type == self._LP.TYPE_CORINA or to_type == self._LP.TYPE_LIGPREP):
                    buffer.append(Ligand(smile=mol.get_smile(), ligand_number=mol.get_ligand_number(),
                                         original_smile=mol.get_original_smile(),
                                         enumeration=mol.get_enumeration(),
                                         molecule=OpenEyeMolToRDkitMol(mol.get_molecule(), bySMILES=bySMILES),
                                         mol_type=self._LP.TYPE_RDKIT,
                                         name=mol.get_name()))
                elif (from_type in [self._LP.TYPE_RDKIT, self._LP.TYPE_LIGPREP, self._LP.TYPE_CORINA, self._LP.TYPE_OMEGA]) and \
                (to_type in [self._LP.TYPE_RDKIT, self._LP.TYPE_LIGPREP, self._LP.TYPE_CORINA, self._LP.TYPE_OMEGA]):
                    buffer.append(Ligand(smile=mol.get_smile(), ligand_number=mol.get_ligand_number(),
                                         original_smile=mol.get_original_smile(),
                                         enumeration=mol.get_enumeration(),
                                         molecule=mol.get_molecule(),
                                         mol_type=to_type,
                                         name=mol.get_name()))
                else:
                    raise ValueError(f"The translation of {from_type} to {to_type} is not yet supported.")
            return buffer
