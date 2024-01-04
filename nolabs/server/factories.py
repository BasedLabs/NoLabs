def api_handlers_factory(is_test=False, is_demo=False):
    # assert (not is_test and is_demo) or (is_test and not is_demo)
    from nolabs.server.api_handlers.conformations import ConformationsApiHandler
    from nolabs.server.api_handlers.protein_design import ProteinDesignApiHandler

    if is_demo:
        from nolabs.server.api_handlers.amino_acid_demo import AminoAcidLabDemoApiHandler
        from nolabs.server.api_handlers.drug_target_demo import DrugTargetDemoApiHandler

        return AminoAcidLabDemoApiHandler(), DrugTargetDemoApiHandler(), ConformationsApiHandler(), ProteinDesignApiHandler()

    # if is_test:
    #    return AminoAcidLabApiMockHandler(), DrugTargetApiMockHandler()

    from nolabs.server.api_handlers.amino_acid import AminoAcidLabApiHandler
    from nolabs.server.api_handlers.drug_target import DrugTargetApiHandler
    return AminoAcidLabApiHandler(), DrugTargetApiHandler(), ConformationsApiHandler(), ProteinDesignApiHandler()
