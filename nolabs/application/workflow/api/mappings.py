from typing import Union

from nolabs.application.workflow.api.schema import ItemsSchema, PropertySchema
from nolabs.application.workflow.component import Items, Parameter, Property


def map_items(i: Items, schema: Parameter) -> ItemsSchema:
    ref = None

    if i.ref:
        ref = i.ref

    if i.items and i.items.ref:
        ref = i.items.ref

    definition = None

    if ref:
        definition = schema.defs[schema.get_ref_type_name(ref)]

    items = None

    if i.items:
        items = (
            map_items(i.items, schema)
            if not isinstance(i.items, list)
            else [map_items(item, schema) for item in i.items]
        )

    properties = (
        {
            name: map_property(prop, schema)
            for name, prop in (
                definition.properties if definition.properties else {}
            ).items()
        }
        if definition
        else {
            name: map_property(prop, schema)
            for name, prop in (i.properties if i.properties else {}).items()
        }
    )

    return ItemsSchema(
        type=i.type,
        properties=properties if not items else None,
        required=i.required,
        description=i.description,
        enum=i.enum,
        const=i.const,
        format=i.format,
        default=i.default,
        example=i.example,
        items=items,
    )


def map_property(p: Union[Property, Items], schema: Parameter) -> PropertySchema:
    ref = None

    if p.ref:
        ref = p.ref

    if p.items and p.items.ref:
        ref = p.items.ref

    definition = None

    if ref:
        definition = schema.defs[schema.get_ref_type_name(ref)]

    items = None

    if p.items:
        items = (
            map_items(p.items, schema)
            if not isinstance(p.items, list)
            else [map_items(item, schema) for item in p.items]
        )

    properties = (
        {
            name: map_property(prop, schema)
            for name, prop in (
                definition.properties if definition.properties else {}
            ).items()
        }
        if definition
        else {
            name: map_property(prop, schema)
            for name, prop in (p.properties if p.properties else {}).items()
        }
    )

    return PropertySchema(
        type=p.type,
        properties=properties if not items else None,
        required=p.required,
        description=p.description,
        enum=p.enum,
        const=p.const,
        format=p.format,
        default=p.default,
        example=p.example,
        title=p.title,
        anyOf=[map_property(prop, schema=schema) for prop in p.anyOf],
        ref=ref,
        items=items,
    )
