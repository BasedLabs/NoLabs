from copy import deepcopy
from typing import Optional, List

from typing_extensions import Literal

from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

from dockstream.core.stereo_enumerator import StereoEnumerator
from dockstream.core.ligand.ligand import Ligand, get_next_enumeration_number_for_ligand


class RDKitStereoEnumeratorParameters(BaseModel):
    try_embedding: bool = True
    unique: bool = True
    max_isomers: int = 1024
    rand: Optional[int] = 0xf00d


class RDKitStereoEnumerator(StereoEnumerator, BaseModel):

    backend: Literal["RDKit"] = "RDKit"
    parameters: RDKitStereoEnumeratorParameters = RDKitStereoEnumeratorParameters()

    def __init__(self, **data):
        super().__init__(**data)

    def enumerate(self, ligands: List[Ligand]) -> List[Ligand]:
        new_ligands_list = []
        opts = StereoEnumerationOptions(tryEmbedding=self.parameters.try_embedding,
                                        unique=self.parameters.unique,
                                        maxIsomers=self.parameters.max_isomers,
                                        rand=self.parameters.rand)
        for ligand in ligands:
            molecule = Chem.MolFromSmiles(ligand.get_smile())
            if not molecule:
                # could not build molecule, keep the original ligand
                new_ligands_list.append(deepcopy(ligand))
                continue

            isomers = tuple(EnumerateStereoisomers(molecule, options=opts))
            if len(isomers) == 0:
                # could not enumerate, keep original ligand
                new_ligands_list.append(deepcopy(ligand))
                continue

            # loop over stereo-isomers, translate them into smiles and create new ligand objects from them
            for new_smile_id, new_smile in enumerate(sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers)):
                new_ligands_list.append(Ligand(smile=new_smile,
                                               original_smile=ligand.get_original_smile(),
                                               ligand_number=ligand.get_ligand_number(),
                                               enumeration=get_next_enumeration_number_for_ligand(new_ligands_list,
                                                                                                  ligand.get_ligand_number()),
                                               molecule=None,
                                               mol_type=None,
                                               name=ligand.get_name()))
        return new_ligands_list
