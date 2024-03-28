
class TagAdditionsEnum:

    TAG_NAME = "name"
    TAG_SMILES = "smiles"
    TAG_ORIGINAL_SMILES = "original_smiles"
    TAG_LIGAND_ID = "ligand_id"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")