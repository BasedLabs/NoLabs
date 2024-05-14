from typing import Type, Dict, Any

from pydantic import BaseModel


class Mapper:
    def __init__(self, target_model: Type[BaseModel]):
        self.target_model = target_model
        self.mappings: Dict[str, str] = {}

    def add_mapping(self, target_field: str, source_model: Type[BaseModel], source_field: str):
        self.mappings[target_field] = (source_model, source_field)

    def validate_and_map(self, *source_objects: BaseModel) -> BaseModel:
        target_data = {}
        errors = []

        for name, field in self.target_model.__fields__:
            if name in self.mappings:
                source_model, source_field = self.mappings[name]
                source_obj = next((obj for obj in source_objects if isinstance(obj, source_model)), None)

                if source_obj is None:
                    errors.append(f"Missing source object for field '{name}'")
                    continue

                if not hasattr(source_obj, source_field):
                    errors.append(f"Missing field '{source_field}' in source object '{source_model.__name__}'")
                    continue

                value = getattr(source_obj, source_field)
                if not isinstance(value, field.annotation):
                    errors.append(f"Type mismatch for field '{name}': expected {field.annotation}, got {type(value)}")
                    continue

                target_data[name] = value
            else:
                errors.append(f"Missing mapping for target field '{name}'")

        if errors:
            raise ValueError("Mapping validation errors: " + ", ".join(errors))

        return self.target_model(**target_data)


from pydantic import BaseModel

class Source1(BaseModel):
    a: int
    b: str

class Source2(BaseModel):
    c: float

class Target(BaseModel):
    res: int
    res2: float


# Define source objects
source1 = Source1(a=5, b="example")
source2 = Source2(c=3.14)

# Set up the mapper
mapper = Mapper(Target)
mapper.add_mapping('res', Source1, 'a')
mapper.add_mapping('res2', Source2, 'c')

# Perform validation and mapping
try:
    target_instance = mapper.validate_and_map(source1, source2)
    print(target_instance)
except ValueError as e:
    print(e)