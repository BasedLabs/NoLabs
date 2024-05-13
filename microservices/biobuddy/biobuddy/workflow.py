from typing import List


class Component:
    def __init__(self, name: str, inputs: List[str] = None, outputs: List[str] = None, description: str = None):
        self.name = name
        self.inputs = inputs or []
        self.outputs = outputs or []
        self.description = description

# Creating 10 fake components related to biology

components = [
    Component(
        name="DownloadFromRCSB",
        outputs=["fasta"],
        description="Downloads data from RCSB PDB in FASTA format."
    ),
    Component(
        name="DownloadFromChembl",
        outputs=["smiles"],
        description="Downloads ligands from ChEMBL in SMILES format."
    ),
    Component(
        name="DockingProteinOnLigand",
        inputs=["pdb", "smiles"],
        outputs=["pdb"],
        description="Performs docking of protein on ligand, accepts PDB file and SMILES string, and generates a PDB file."
    ),
    Component(
        name="DockingProteinOnProtein",
        inputs=["pdb1", "pdb2"],
        outputs=["docked_pdb"],
        description="Performs docking of one protein on another, accepts two PDB files, and generates one PDB file."
    ),
    Component(
        name="PredictSolubility",
        inputs=["protein_sequence"],
        outputs=["solubility_score"],
        description="Predicts solubility of a protein, accepts protein sequence, and returns a float solubility score."
    ),
    Component(
        name="CalculateBlast",
        inputs=["sequence"],
        outputs=["blast_results"],
        description="Performs BLAST (Basic Local Alignment Search Tool), accepts a sequence, and returns BLAST results."
    ),
    Component(
        name="GenerateProteins",
        outputs=["pdb"],
        description="Generates proteins based on input conditions and returns a PDB file."
    ),
    Component(
        name="RunFolding",
        inputs=["protein_sequence"],
        outputs=["folded_pdb"],
        description="Runs folding of proteins, accepts protein sequence, and returns a folded PDB file."
    ),
    # Add more components here
]



