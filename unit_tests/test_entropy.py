import sys
import unittest
import os
sys.path.insert(1, '/home/radeksz/bioinf/uni-gem-rec/MSAquality')
from unittest import TestCase
from entropy import Entropy

class TestEntropy(TestCase):

    def testReturnEntropy(self):
        cwd = os.getcwd()
        entropy = Entropy(cwd, "testSave", 0.8)
        testDict = {"id1": "MMMSSS"}


if __name__ == "__main__":
    unittest.main()