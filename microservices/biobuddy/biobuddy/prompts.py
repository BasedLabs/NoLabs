from biobuddy.workflow import components


def generate_system_prompt() -> str:
    return ("You are a biobuddy. You are a research assistant that helps creating new drugs, researching biology, "
        "pulling information from different sources and running experiments. "
        "If a user asks you to pull the targets only do so for proteins and use rcsb pdb. If a user asks for "
        "drugs/ligands then use chembl.")


def generate_strategy_prompt(tools_description: str, query: str) -> str:
    return (f"Given the available tools ({tools_description}), decide whether a direct reply is sufficient "
        f"or if a specific plan involving these tools is necessary. "
            f"If a calling function calls are not needed, provide the direct reply (without stating that it's a direct reply, just the text of the reply). "
        f"If function calls are needed, reply with the list (with a tag <PLAN> before the list) of queries for function calls in the format: \"<PLAN>: ['Call function_1 in order to do X', 'Call function_2 to achieve Y']\". "
            f"If you are asked to generate a workflow just reply '<WORKFLOW>' tag. I'll process the query separately."
        f"If you are asked SPECIFICALLY to do the literature research on some topic, write only tag <RESEARCH>. "
        f"\n\nQuery: \"{query}\"\n\n")

workflow_json = '''
{
  "workflow": {
    "name": "My Biology Workflow",
    "nodes": [
      {"id": "1", "name": "DownloadFromRCSB", "type": "component", "inputs": [], "outputs": ["fasta_file"], "description": "Downloads data from RCSB PDB in FASTA format."},
      {"id": "2", "name": "DownloadFromChembl", "type": "component", "inputs": [], "outputs": ["smiles_string"], "description": "Downloads ligands from ChEMBL in SMILES format."},
      {"id": "3", "name": "DockingProteinOnLigand", "type": "component", "inputs": ["pdb_file", "smiles_string"], "outputs": ["pdb_file"], "description": "Performs docking of protein on ligand, accepts PDB file and SMILES string, and generates a PDB file."},
      {"id": "4", "name": "PredictSolubility", "type": "component", "inputs": ["fasta_file"], "outputs": ["solubility_score"], "description": "Predicts solubility of a protein, accepts protein sequence, and returns a float solubility score."},
      {"id": "5", "name": "RunFolding", "type": "component", "inputs": ["fasta_file"], "outputs": ["pdb_file"], "description": "Runs folding of proteins, accepts protein sequence, and returns a folded PDB file."}
    ],
    "edges": [
      {"from": "1", "to": "3"},
      {"from": "2", "to": "3"},
      {"from": "1", "to": "4"},
      {"from": "1", "to": "5"}
    ]
  }
}
'''
def generate_workflow_prompt(query: str) -> str:
    return (f"Generate a JSON which would be used to construct a workflow graph. Reply just by generating the json, don't add any explanations or comments. "
            f"Workflow consists of components. A hypothetical example of the output could be:"
            f"{workflow_json}"
            f"Available components to choose from are: {', '.join(component.name for component in components)}"
            f"\n\nUser Query: \"{query}\"\n\n")