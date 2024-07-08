from biobuddy.workflow import components


def generate_system_prompt() -> str:
    return ("You are a biobuddy. You are a research assistant that helps creating new drugs, researching biology, "
        "pulling information from different sources and running experiments. "
        "If a user asks you to pull the targets only do so for proteins and use rcsb pdb. If a user asks for "
        "drugs/ligands then use chembl.")


def generate_strategy_prompt(tools_description: str, query: str, available_components: str, current_workflow: str) -> str:
    return (f"Given the available tools ({tools_description}), decide whether a direct reply is sufficient "
            f"or if a specific plan involving these tools is necessary. "
            f"If calling function calls are not needed, provide the direct reply (without stating that it's a direct reply, just the text of the reply). "
            f"If the query explicitly asks to pull specific data (such as named proteins or ligands), reply with a plan of all actions required based on the query and tools descriptions in the format (have as many actions as needed): "
            f"\"<ACTION> 1. Call tool_1 in order to do X <END_ACTION> <ACTION> 2. Call tool_2 to achieve Y <END_ACTION> <ACTION> 3. Call tool_3 to do something else <END_ACTION>\". "
            f"Be careful to actually use all the functions, i.e. if a user asks to pull two specific molecules and it can't be done with just one action, do it with two actions. "
            f"DON'T GENERATE THE PLAN IF IT'S NOT NEEDED, i.e. if a user asks to just adjust a workflow or has a simple general question. "
            f"DO NOT GENERATE A PLAN WITH ACTIONS IF THE USER DOES NOT SPECIFY THE PROTEINS OR LIGANDS THEY WANT, such as in a request to 'generate a workflow with some proteins and ligands.' "
            f"If you are asked SPECIFICALLY to do the literature research on some topic, write only tag <RESEARCH>. "
            f"\n\nAdditionally, if the workflow is empty or needs adjustments based on the query, add or adjust components from the available components. "
            f"Ensure to connect the components appropriately using the connections format provided. If new components are added, assign them ids like newId_1, newId_2, etc., and encapsulate the entire workflow in <WORKFLOW> and <END_WORKFLOW> tags.\n"
            f"For instance, if you are asked to fold a protein and there's no workflow component available for folding, you should connect proteins holder to esmfold light by default, if a particular method is not specified."
            f"\n\nQuery: \"{query}\"\n\n"
            f"Available Components: {available_components}\n\n"
            f"Current Workflow: {current_workflow}\n\n"
            f"Like if the query tells to do something with proteins or ligands, make sure that these data holder components are present."
            f"Your response should include the plan and the adjusted workflow if necessary. Add a short summary of changes before the JSON if you are adjusting the workflow (but don't add 'Adjusted workflow:' before the json. The workflow should be in the JSON format: "
            f"<WORKFLOW> {{\"workflow_components\": [{{\"id\": \"newId_1\", \"name\": \"component_name\", \"description\": \"component_description\", \"connections\": [{{\"source_output\": [\"source_output_name\"], \"target_input\": [\"input_name\"], \"source_component_id\": \"source_id\"}}]}}]}} <END_WORKFLOW>."
            f"input_name and output_name come from the list of inputs and outputs specified in the component. For instance, esmfold_light has 'proteins_with_fasta' as input and 'proteins_with_pdb' as output. So source_component_id should be the id of proteins component,"
            f" the source_output should be 'proteins' since it's an output for proteins component and target_input should be 'proteins_with_fasta' "
            f"Example of a good workflow freshly generated for docking: "
            f"{{\"workflow_components\": ["
            f"{{\"id\": \"newId_1\", \"name\": \"Proteins\", \"description\": \"Proteins datasource\", \"connections\": []}}, "
            f"{{\"id\": \"newId_2\", \"name\": \"Esmfold light\", \"description\": \"Protein folding using Esmfold light\", \"connections\": ["
            f"{{\"source_output\": [\"proteins\"], \"target_input\": [\"proteins_with_fasta\"], \"source_component_id\": \"newId_1\"}}]}}, "
            f"{{\"id\": \"newId_3\", \"name\": \"Ligands\", \"description\": \"Ligands datasource\", \"connections\": []}}, "
            f"{{\"id\": \"newId_4\", \"name\": \"DiffDock\", \"description\": \"Prediction of protein-ligand complexes\", \"connections\": ["
            f"{{\"source_output\": [\"proteins_with_pdb\"], \"target_input\": [\"proteins_with_pdb\"], \"source_component_id\": \"newId_2\"}}, "
            f"{{\"source_output\": [\"ligands\"], \"target_input\": [\"ligands\"], \"source_component_id\": \"newId_3\"}}]}}]}}")


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