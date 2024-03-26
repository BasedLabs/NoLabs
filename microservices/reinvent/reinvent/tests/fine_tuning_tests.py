import os
import unittest

from microservice import FineTuning


class FineTuningTests(unittest.TestCase):
    path = os.path.dirname(os.path.abspath(__file__))

    def get_fine_tuning(self) -> FineTuning:
        return FineTuning(os.path.join(self.path, 'finetuning'))

    def test_upper(self):
        fine_tuning = FineTuning(self.path)



if __name__ == '__main__':
    unittest.main()
