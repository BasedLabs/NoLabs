import json
from collections import OrderedDict
from typing import Dict, Any

from dockstream.core.input_model import AzdockInput
from dockstream.utils.schema import (
    remove_schema_properties,
    add_boolean_guards_for_schema_properties,
    replacekey,
    addsibling,
    delsibling,
    copytitle,
    makeconst,
)


def patch_schema_generic(schema):

    # Replace "anyOf" with "oneOf".
    schema = replacekey(schema)

    # Add "type": "object" to any elements that contain "oneOf": [...].
    schema = addsibling(schema)

    # Delete "type": "string" for "enum".
    schema = delsibling(schema, {"enum": "type"})

    # Delete most of the stuff for "const".
    schema = delsibling(schema, {"const": "type"})
    schema = delsibling(schema, {"const": "default"})
    schema = delsibling(schema, {"const": "title"})

    # Copy title from $refs into oneOf.
    schema = copytitle(schema, schema)

    return schema


def add_aws_scp_choice_for_schrodinger(schema, aws, scp) -> None:
    """Add AWS/SCP choice."""
    props = schema.get("properties", {})
    props["use_aws_or_scp"] = {
        "title": "Execute on AWS or SCP",
        "default": "AWS",
        "enum": ["AWS", "SCP",],
    }
    if "dependencies" not in schema:
        schema["dependencies"] = {}

    aws_props = makeconst(aws)
    aws_props["properties"]["use_aws_or_scp"] = {"enum": ["AWS"]}
    aws_props.pop("type", None)

    scp_props = makeconst(scp)
    scp_props["properties"]["use_aws_or_scp"] = {"enum": ["SCP"]}
    scp_props.pop("type", None)

    schema["dependencies"]["use_aws_or_scp"] = {
        "oneOf": [aws_props, scp_props],
        "type": "object",
    }


def schema_extra_glideparams(topschema: Dict[str, Any]) -> None:
    schema = topschema.get("definitions", {}).get("GlideParameters", {})
    remove_schema_properties(
        schema,
        [
            "prefix_execution",
            "binary_location",
            "time_limit_per_compound",
            "token_guard",
            # "parallelization",
            "glide_flags",
        ],
    )

    aws = {
        "prefix_execution": "module load schrodinger/2020-4-js-aws",
        "glide_flags": {"-HOST": "cpu-only"},
    }

    scp = {
        "prefix_execution": "module load schrodinger/2019-4",
        "glide_flags": {"-HOST": "localhost"},
        "token_guard": {
            "prefix_execution": "module load schrodinger/2019-4",
            "token_pools": {"GLIDE_SP_DOCKING": 32, "GLIDE_SUITE_11DEC2020": 32},
            "wait_interval_seconds": 30,
            "wait_limit_seconds": 900,
        },
    }

    add_aws_scp_choice_for_schrodinger(schema, aws, scp)

    add_boolean_guards_for_schema_properties(schema, ["advanced_glide_keywords"])

    schema.pop("title", None)
    schema.pop("description", None)


def schema_extra_advanced_glide_keywords(topschema):
    schema = topschema.get("definitions", {}).get("AdvancedGlideKeywords", {})
    props = schema.get("properties", {})
    props.get("maestro_file", {})["format"] = "uri"
    props.get("REF_LIGAND_FILE", {})["format"] = "uri"


def schema_extra_glidekw(topschema):
    schema = topschema.get("definitions", {}).get("GlideKeywords", {})

    # Glide keywords override Advanced Glide Keywords > Maestro File.
    # We show all keywords, otherwise their default values
    # will always override whatever is in the Maestro file,
    # without possibility to change those values from the GUI.
    remove_schema_properties(
        schema,
        [
            # "AMIDE_MODE",
            # "POSTDOCK_NPOSE",
            # "POSTDOCKSTRAIN",
            "LIGANDFILE",  # This keyword is overridden by DockStream input.
            "POSE_OUTTYPE",  # This keyword is overridden by DockStream (single outtype).
        ],
    )

    # Change gridfile "format" to "uri".
    props = schema.get("properties", {})
    gridfile = props.get("GRIDFILE", {})
    gridfile.pop("oneOf", {})  # Remove possibility for multiple gridfiles.
    gridfile.pop("anyOf", {})  # Remove possibility for multiple gridfiles.
    gridfile["type"] = "string"
    gridfile["format"] = "uri"


def schema_extra_ligprepparams(topschema):
    schema = topschema.get("definitions", {}).get("LigprepLigandPreparatorParameters", {})
    remove_schema_properties(
        schema,
        [
            "prefix_execution",
            "binary_location",
            "token_guard",
            # "parallelization",
            "filter_file",
            "command_line_parameters",
        ],
    )

    aws = {
        "prefix_execution": "module load schrodinger/2020-4-js-aws",
        "command_line_parameters": {"-HOST": "cpu-only",},
    }

    scp = {
        "prefix_execution": "module load schrodinger/2019-4",
        "command_line_parameters": {"-HOST": "localhost"},
        "token_guard": {
            "prefix_execution": "module load schrodinger/2019-4",
            "token_pools": {"GLIDE_SP_DOCKING": 32, "GLIDE_SUITE_11DEC2020": 32},
            "wait_interval_seconds": 30,
            "wait_limit_seconds": 900,
        },
    }

    add_aws_scp_choice_for_schrodinger(schema, aws, scp)

    add_boolean_guards_for_schema_properties(schema, ["use_epik", "chirality"])

    # Switch default for IGNORE_CHIRALITIES to True.
    # IGNORE_CHIRALITIES was added for the GUI,
    # users get this functionality from DockStream by passing `{"-ac": ""}`
    # as an element in command_line_arguments.
    # If IGNORE_CHIRALITIES is absent from the model,
    # it should be False (do not add "-ac"),
    # thus, for backward compatibility,
    # the default value for it in the Field specification is False.
    # However, default value in the GUI should be set to True,
    # and we change to True here in the schema.
    schema.get("properties", {}).get("IGNORE_CHIRALITIES", {})["default"] = True


def schema_extra_DockingInput(topschema):
    schema = topschema.get("definitions", {}).get("DockingInput", {})
    remove_schema_properties(schema, ["header",])
    schema.pop("title", None)
    schema.pop("description", None)


def schema_extra_ligand_preparator(topschema):
    schema = topschema.get("definitions", {}).get("LigprepLigandPreparator", {})
    remove_schema_properties(
        schema, ["align", "output", "ligands"]  # "pool_id", "input"
    )
    pool_id = schema.get("properties", {}).get("pool_id", {})
    pool_id.pop("type", {})
    pool_id["const"] = "pool"

    schema.get("properties", {}).get("parameters", {}).pop("default", {})


def schema_extra_ligprep_input(topschema):
    schema = topschema.get("definitions", {}).get("Input", {})
    schema["properties"] = {"type": {"const": "console"}}
    # Previous assignment deletes all other properties.

    schema.pop("title", None)
    schema.pop("description", None)


def schema_extra_docker(topschema):
    schema = topschema.get("definitions", {}).get("Glide", {})
    remove_schema_properties(schema, ["ligands"])  # "input_pools"

    add_boolean_guards_for_schema_properties(schema, ["output"])

    input_pools = schema.get("properties", {}).get("input_pools", {})
    input_pools.pop("type", {})
    input_pools["const"] = "pool"

    run_id = schema.get("properties", {}).get("run_id", {})
    run_id.pop("type", {})
    run_id["const"] = "docking_run"


def schema_extra_glide(topschema):
    schema = topschema.get("definitions", {}).get("Glide", {})
    schema.get("properties", {}).get("parameters", {}).pop("default", None)


def schema_extra_poses(topschema):
    schema = topschema.get("definitions", {}).get("Poses", {})
    remove_schema_properties(schema, ["overwrite"])  # "poses_path"

    path = schema.get("properties", {}).get("poses_path", {})
    path.pop("type", {})
    path["const"] = "{{component.path}}/output/poses.sdf"


def schema_extra_scores(topschema):
    schema = topschema.get("definitions", {}).get("Scores", {})
    remove_schema_properties(schema, ["overwrite"])  # "scores_path"

    path = schema.get("properties", {}).get("scores_path", {})
    path.pop("type", {})
    path["const"] = "{{component.path}}/output/scores.csv"


def schema_extra_LigandPreparation(topschema: Dict):
    schema = topschema.get("definitions", {}).get("LigandPreparation", {})
    schema.pop("title", None)
    schema.pop("description", None)


def schema_extra_AzdockInput(schema: Dict):
    schema["title"] = "Docking"
    schema.pop("description", None)
    ordered_props = OrderedDict(schema["properties"])
    ordered_props.update(
        {
            "name": {"type": "string", "title": "Docking name"},
            "description": {"type": "string", "title": "Docking description"},
        }
    )
    ordered_props.move_to_end("description", last=False)
    ordered_props.move_to_end("name", last=False)
    schema["properties"] = ordered_props
    schema["required"].append("name")


def patch_schema_dockstream(schema):

    schema_extra_ligprepparams(schema)

    schema_extra_ligand_preparator(schema)

    schema_extra_ligprep_input(schema)

    schema_extra_glidekw(schema)

    schema_extra_advanced_glide_keywords(schema)

    schema_extra_glideparams(schema)

    # Glide class inherits fields from Docker.
    schema_extra_docker(schema)
    schema_extra_glide(schema)

    schema_extra_scores(schema)

    schema_extra_poses(schema)

    schema_extra_DockingInput(schema)

    schema_extra_LigandPreparation(schema)

    schema_extra_AzdockInput(schema)  # root

    return schema


def main():
    schema = AzdockInput.schema()
    schema = patch_schema_dockstream(schema)
    schema = patch_schema_generic(schema)
    print(json.dumps(schema, indent=2))


if __name__ == "__main__":
    main()
