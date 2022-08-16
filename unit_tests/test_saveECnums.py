import unittest
import sys
import os
sys.path.insert(1, '/home/radeksz/bioinf/uni-gem-rec/ecnumModel')
from unittest import TestCase
import saveECnums


class SaveEcNumsTest(TestCase):

    def testChangePathToSave(self):
        ecNumSaver = saveECnums.ECnumSaver()
        path = ecNumSaver.path
        ecNumSaver.changePathToSave("newPah")
        self.assertNotEqual(path, ecNumSaver.path)

    def testSavingMethod(self):
        ecNumSaver = saveECnums.ECnumSaver()
        path = os.path.join(ecNumSaver.path, "testKoEcnum.fasta")
        self.assertFalse(os.path.isfile(path))
        ecNumSaver.saveToFile("testSeq", "testKoEcnum")
        self.assertTrue(os.path.isfile(path))
        os.remove(path)



if __name__ == "__main__":
    unittest.main()