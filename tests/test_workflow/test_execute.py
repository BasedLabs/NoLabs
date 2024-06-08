from unittest import IsolatedAsyncioTestCase

from tests.tests_preparations import mongo_connect


class TestWorkflowExecute(IsolatedAsyncioTestCase):
    def setUp(self):
        mongo_connect()

    def test_linear_execute(self):
        pass