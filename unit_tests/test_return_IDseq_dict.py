import unittest
import sys
import os
sys.path.insert(1, '/home/radeksz/bioinf/uni-gem-rec/ecnumModel')
from unittest import TestCase
import return_IDseq_dict

class SaveEcNumsTest(TestCase):

    def testReturnigTypeAndSavingToFile(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFasta:
            testFasta.write(">test:testID")
            testFasta.write("\n")
            testFasta.write("testSeq")
        with open(testPath) as testFasta:
            allLines = testFasta.readlines()
            print(allLines)
            self.assertNotEqual(len(allLines), 0)
        returnIdSeq_dict = return_IDseq_dict.ReturnIndexSequence(testPath)
        idSeq_dict = returnIdSeq_dict.returnPeptideDict()
        dictBool = False
        if dict == type(idSeq_dict):
            dictBool = True
        self.assertTrue(dictBool)
        os.remove(testPath)
    
if __name__ == "__main__":
    unittest.main()