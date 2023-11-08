
def api_handlers_factory(is_test=False, is_demo=False):
    #assert (not is_test and is_demo) or (is_test and not is_demo)

    if is_demo:
        from src.server.api_handlers.amino_acid_demo import AminoAcidLabDemoApiHandler
        from src.server.api_handlers.drug_target_demo import DrugTargetDemoApiHandler

        return AminoAcidLabDemoApiHandler(), DrugTargetDemoApiHandler()

    # if is_test:
    #    return AminoAcidLabApiMockHandler(), DrugTargetApiMockHandler()

    from src.server.api_handlers.amino_acid import AminoAcidLabApiHandler
    from src.server.api_handlers.drug_target import DrugTargetApiHandler
    return AminoAcidLabApiHandler(), DrugTargetApiHandler()
