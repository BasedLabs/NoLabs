from dockstream.containers.container import ConfigurationContainer


class AnalysisContainer(ConfigurationContainer):
    """Class that takes a lot of arguments in the form of a JSON configuration (as dictionary, string or file path)
       and, optionally, performs a JSON Schema validation."""

    def __init__(self, conf, validation=True):
        super().__init__(conf=conf)

        # TODO: include validation with JSON Schema
        if validation:
            self.validate()

    def validate(self):
        pass
        # load the Schema
        #path = os.path.join(files_paths.move_up_directory(__file__, 1),
        #                    "docking", "json_schemas",
        #                    "BuildingConfiguration.json")
        #schema = load_json.loadJSON(path=path)

        # instantiate validator and perform the check
        #validator = Draft7Validator(schema=schema)
        #validator.validate(self._conf)
