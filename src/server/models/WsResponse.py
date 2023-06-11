class WsResponseBase:
    def to_json(self):
        pass


class WsProgressResponse(WsResponseBase):
    def __init__(self, progress):
        self.progress = progress

    def to_json(self):
        return {
            'progress': self.progress
        }


class WsResultResponse(WsResponseBase):
    def __init__(self, sequence, localisation, folding):
        self.sequence = sequence
        self.localisation = localisation
        self.folding = folding

    def to_json(self):
        return {
            'sequence': self.sequence,
            'localisation': self.localisation,
            'folding': self.folding
        }
