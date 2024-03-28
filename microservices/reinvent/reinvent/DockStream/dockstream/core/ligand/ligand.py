from copy import deepcopy
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.tag_additions_enum import TagAdditionsEnum


class Ligand:
    """This class bundles all information on a ligand, including all molecule instances present."""

    def __init__(self, smile: str, ligand_number: int, enumeration=0, molecule=None, mol_type=None, name=None, original_smile=None):
        # initialize
        self._LP = LigandPreparationEnum()
        self._TA = TagAdditionsEnum()
        self._known_types = [self._LP.TYPE_RDKIT, self._LP.TYPE_OPENEYE, self._LP.TYPE_OMEGA,
                             self._LP.TYPE_CORINA, self._LP.TYPE_GOLD, self._LP.TYPE_LIGPREP, None]
        # set attributes
        self._smile = self._check_smile(smile)
        self._original_smile = original_smile
        self._ligand_number = self._check_ligand_number(ligand_number)
        self._enumeration = self._check_enumeration(enumeration)
        self._molecule = molecule
        self._mol_type = self._check_mol_type(mol_type)
        self._name = name
        self._conformers = []

    def __repr__(self):
        return "<Ligand id: %s, enumeration: %s, smile: %s>" % (self.get_ligand_number(), self.get_enumeration(), self.get_smile())

    def __str__(self):
        return f"Ligand id: {self.get_ligand_number()}, enumeration: {self.get_enumeration()}, " + \
               f"name: {self.get_name()}, smile: {self.get_smile()}, original_smile: {self.get_original_smile()}, " + \
               f"mol_type: {self.get_mol_type()}, has molecule: {True if self.get_molecule() is not None else False}."

    def get_clone(self):
        clone = Ligand(smile=self.get_smile(),
                       ligand_number=self.get_ligand_number(),
                       enumeration=self.get_enumeration(),
                       molecule=deepcopy(self.get_molecule()),
                       mol_type=self.get_mol_type(),
                       name=self.get_name(),
                       original_smile=self.get_original_smile())
        for conformer in self.get_conformers():
            clone.add_conformer(deepcopy(conformer))
        return clone

    def __copy__(self):
        return self.get_clone()

    def __deepcopy__(self, memo):
        return self.get_clone()

    def set_name(self, name: str):
        self._name = name

    def get_name(self) -> str:
        return self._name

    def add_conformer(self, conformer):
        self._conformers.append(conformer)

    def set_conformers(self, conformers: list):
        self._conformers = conformers

    def get_conformers(self):
        return self._conformers

    def clear_conformers(self):
        self._conformers = []

    def set_molecule(self, molecule):
        self._molecule = molecule

    def get_molecule(self):
        return self._molecule

    def _check_mol_type(self, mol_type) -> str:
        if mol_type not in self._known_types:
            raise ValueError(f"Type {mol_type} not in list of supported types.")
        return mol_type

    def set_mol_type(self, mol_type):
        self._mol_type = self._check_mol_type(mol_type)

    def get_mol_type(self):
        return self._mol_type

    def _check_smile(self, smile: str) -> str:
        if not isinstance(smile, str):
            raise ValueError(f"Field smile must be a string not of type {type(smile)}.")
        return smile

    def set_smile(self, smile: str):
        self._smile = self._check_smile(smile)

    def get_smile(self):
        return self._smile

    def set_original_smile(self, smile: str):
        self._original_smile = smile

    def get_original_smile(self):
        return self._original_smile

    def _check_ligand_number(self, ligand_number: int):
        if not isinstance(ligand_number, int) or ligand_number < 0:
            raise ValueError(f"Ligand number must be an integer value (minimally 0), not {ligand_number}.")
        return ligand_number

    def set_ligand_number(self, ligand_number: int):
        self._ligand_number = self._check_ligand_number(ligand_number)

    def get_ligand_number(self):
        return self._ligand_number

    def _check_enumeration(self, enumeration: int):
        if not isinstance(enumeration, int) or enumeration < 0:
            raise ValueError(f"Enumeration must be an integer value (minimally 0), not {enumeration}.")
        return enumeration

    def set_enumeration(self, enumeration: int):
        self._enumeration = self._check_enumeration(enumeration)

    def get_enumeration(self):
        return self._enumeration

    def get_identifier(self):
        return str(self.get_ligand_number()) + ':' + str(self.get_enumeration())

    def _add_title_to_molecule(self, molecule, title):
        if self.get_mol_type() in [self._LP.TYPE_RDKIT, self._LP.TYPE_CORINA, self._LP.TYPE_GOLD, self._LP.TYPE_OMEGA]:
            molecule.SetProp("_Name", str(title))
        elif self.get_mol_type() == self._LP.TYPE_OPENEYE:
            molecule.SetTitle(str(title))

    def _add_tag_to_molecule(self, molecule, tag, value):
        if self.get_mol_type() in [self._LP.TYPE_RDKIT, self._LP.TYPE_CORINA, self._LP.TYPE_GOLD, self._LP.TYPE_LIGPREP,
                                   self._LP.TYPE_OMEGA]:
            molecule.SetProp(tag, str(value))
        elif self.get_mol_type() == self._LP.TYPE_OPENEYE:
            import openeye.oechem as oechem
            oechem.OESetSDData(molecule, tag, str(value))
        else:
            raise ValueError(f"Cannot add tags to conformer type {self.get_mol_type()}.")

    def add_tags_to_conformers(self):
        if len(self.get_conformers()) > 0:
            for conformer_number, conformer in enumerate(self.get_conformers()):
                self._add_title_to_molecule(conformer, self.get_identifier() + ':' + str(conformer_number))
                if self.get_name() is not None:
                    self._add_tag_to_molecule(conformer, self._TA.TAG_NAME, self.get_name())
                self._add_tag_to_molecule(conformer, self._TA.TAG_LIGAND_ID, self.get_ligand_number())
                self._add_tag_to_molecule(conformer, self._TA.TAG_ORIGINAL_SMILES, self.get_original_smile())
                self._add_tag_to_molecule(conformer, self._TA.TAG_SMILES, self.get_smile())

    def add_tags_to_molecule(self):
        if self.get_molecule() is not None:
            self._add_title_to_molecule(self.get_molecule(), self.get_identifier())
            if self.get_name() is not None:
                self._add_tag_to_molecule(self.get_molecule(), self._TA.TAG_NAME, self.get_name())
            self._add_tag_to_molecule(self.get_molecule(), self._TA.TAG_LIGAND_ID, self.get_ligand_number())
            self._add_tag_to_molecule(self.get_molecule(), self._TA.TAG_ORIGINAL_SMILES, self.get_original_smile())
            self._add_tag_to_molecule(self.get_molecule(), self._TA.TAG_SMILES, self.get_smile())


def get_next_enumeration_number_for_ligand(ligands: list, ligand_id: int):
    max_enumeration = -1
    for ligand in ligands:
        if ligand.get_ligand_number() == ligand_id:
            max_enumeration = max(max_enumeration, ligand.get_enumeration())
    return max_enumeration + 1


def get_enumerations_for_ligand(ligands: list, ligand_id: int):
    ligand_enumerations = []
    for ligand in ligands:
        if ligand.get_ligand_number() == ligand_id:
            ligand_enumerations.append(deepcopy(ligand))
    return ligand_enumerations


def reset_enumerations_for_ligands(ligands: list):
    ligand_identifiers = list(set([lig.get_identifier() for lig in ligands]))
    cur_enum_list = {k: 0 for k in ligand_identifiers}
    for lig in ligands:
        cur_id = lig.get_identifier()
        lig.set_enumeration(cur_enum_list[lig.get_identifier()])
        cur_enum_list[cur_id] += 1


def find_ligand(ligands: list, ligand_id: int, enumeration: int = 0):
    for ligand in ligands:
        if ligand.get_ligand_number() == ligand_id and ligand.get_enumeration() == enumeration:
            return ligand
    return None
