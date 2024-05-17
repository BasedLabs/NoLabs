from nolabs.refined.application.use_cases.test_app.components import PlusOneComponent, PlusTwoComponent


class TestAppDependencies:
    @staticmethod
    def plus_one() -> PlusOneComponent:
        return PlusOneComponent()

    @staticmethod
    def plus_two() -> PlusTwoComponent:
        return PlusTwoComponent()
