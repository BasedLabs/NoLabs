

def generate_system_prompt() -> str:
    return ("You are a biobuddy. You are a research assistant that helps creating new drugs, researching biology, "
        "pulling information from different sources and running experiments. "
        "If a user asks you to pull the targets only do so for proteins and use rcsb pdb. If a user asks for "
        "drugs/ligands then use chembl.")


def generate_strategy_prompt(tools_description: str, query: str) -> str:
    return (f"Given the available tools ({tools_description}), decide whether a direct reply is sufficient "
        f"or if a specific plan involving these tools is necessary. If a plan is not needed, provide the direct reply (without stating that it's a direct reply, just the text of the reply). "
        f"If a plan is needed, reply with the list (with a tag <PLAN> before the list) of function calls in the format <PLAN>: ['Call function_1 in order to do X', 'Call function_2 to achieve Y']. "
        f"If you are asked SPECIFICALLY to do the literature research on some topic, write only tag <RESEARCH>. "
        f"\n\nQuery: \"{query}\"\n\n")