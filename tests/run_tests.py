import unittest
import test_workflow

suite = unittest.TestLoader().loadTestsFromModule(test_workflow)
unittest.TextTestRunner().run(suite)