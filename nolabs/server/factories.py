def api_handlers_factory():
    # assert (not is_test and is_demo) or (is_test and not is_demo)
    from nolabs.server.api_handlers.conformations import ConformationsApiHandler
    from nolabs.server.api_handlers.protein_design import ProteinDesignApiHandler
    return ConformationsApiHandler(), ProteinDesignApiHandler()
