import sys
import os
import unittest
sys.path.insert(1, '/home/radeksz/bioinf/uni-gem-rec/ecnumModel/HMMaproach')
from unittest import TestCase
from autoHMMs import HMM

class AutoHmmTest(TestCase):

    def testAddFilesToList(self):
        with open("test.aln", "w") as fileTest:
            pass 
        cwd = os.getcwd()
        autoHMM = HMM(cwd, "testHMMs")
        pathOfAfile = os.path.join(cwd, "test.aln")
        self.assertEqual(len(autoHMM.fileList), 1)
        self.assertEqual(autoHMM.fileList[0], pathOfAfile)
        os.remove("./test.aln")
        os.removedirs(autoHMM.newDirPathToSaveHMMs)
    
    def testCreateNewDirectory(self):
        with open("test.aln", "w") as fileTest:
            pass 
        cwd = os.getcwd()
        autoHMM = HMM(cwd, "testHMMs")
        self.assertTrue(os.path.isdir(autoHMM.newDirPathToSaveHMMs))
        os.remove("./test.aln")
        os.removedirs(autoHMM.newDirPathToSaveHMMs)
    
    def testCreateHMMs(self):
        with open("test[EC:1.1.1.1].aln", "w") as fileTest:
            fileTest.write(""">ECO26_4684
-------------------------AA--------------------FQDA-----A---
-CD------------YII-AFGGGSPIDTAK----------A--VK--------------
---------ILTANP---------AP-STAYS------------------G-------VG
KVTN-PGVPL-VAINTTAGTAAEMTSNAVIIDT---------A-----RQVKE-VI-IDP
NIIPDI-AVDDA---------------------SVML-----------------------
""") 
        cwd = os.getcwd()
        autoHMM = HMM(cwd, "testHMMs")
        self.assertTrue(os.path.isdir(autoHMM.newDirPathToSaveHMMs))
        autoHMM.createHMMs(os.cpu_count())
        pathOfAfile = os.path.join(cwd, "./testHMMs/test[EC:1.1.1.1].hmm")
        self.assertTrue(os.path.isfile(pathOfAfile))
        os.remove("./test[EC:1.1.1.1].aln")
        os.remove("./testHMMs/test[EC:1.1.1.1].hmm")
        os.removedirs(autoHMM.newDirPathToSaveHMMs)

if __name__ == "__main__":
    unittest.main()
