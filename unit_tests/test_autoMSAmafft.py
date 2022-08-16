from re import L
import sys
import unittest
import os
sys.path.insert(1, '/home/radeksz/bioinf/uni-gem-rec/ecnumModel/HMMaproach')
from unittest import TestCase
from autoMSAmafft import MSA

class TestAutoMSAmafft(TestCase):
    
    def testAddFilesToList(self):
        with open("test.fasta", "w") as fileTest:
            pass 
        cwd = os.getcwd()
        autoMSAmafft = MSA(cwd, "testMSAs")
        pathOfAfile = os.path.join(cwd, "test.fasta")
        self.assertEqual(len(autoMSAmafft.fileList), 1)
        self.assertEqual(autoMSAmafft.fileList[0], pathOfAfile)
        os.remove("./test.fasta")
        os.removedirs(autoMSAmafft.newDirectoryName)
    
    def testCreateNewDirectory(self):
        with open("test.fasta", "w") as fileTest:
            pass 
        cwd = os.getcwd()
        autoMSAmafft = MSA(cwd, "testMSAs")
        self.assertTrue(os.path.isdir(autoMSAmafft.newDirectoryName))
        os.remove("./test.fasta")
        os.removedirs(autoMSAmafft.newDirectoryName)
    

    def testCreateMSAs(self):
        with open("test[EC:1.1.1.1].fasta", "w") as fileTest:
            fileTest.write(""">ECO26_4684
MSFMLALPKISLHGAGAIGDMVNLVAHKQWGKALIVTDGQLVKLGLLDSLFIALDQYEMSYHLFDEVFPNPTEELVQKGYAAF
QDAACDYIIAFGGGSPIDTAKAVKILTANPAPSTAYSGVGKVTNPGVPLVAINTTAGTAAEMTSNAVIIDTARQVKEVIIDPN
IIPDIAVDDASVMLEIPASVTAATGMDALTHAVEAYVSVGAHPLTDANALEAIRLITLWLPKAVDDGHNLEAREQMAFAQYLA
GMAFNSAGLGLVHALAHQPGATHNLPHGVCNAILLPVIENFNRPNAIARFARIAQAMGVDTRDMSEEAASTEAINAIRALSHR
VGIPAGFSQLGITKEDIEGWLDKALADPCVPCNPRTASRDEVRELYLEAL
""") 
        cwd = os.getcwd()
        autoMSAmafft = MSA(cwd, "testHMMs")
        self.assertTrue(os.path.isdir(autoMSAmafft.newDirectoryName))
        autoMSAmafft.createMSAs(os.cpu_count(), 62)
        pathOfAfile = os.path.join(cwd, "./testHMMs/test[EC:1.1.1.1].aln")
        print("#"*80)
        print(pathOfAfile)
        self.assertTrue(os.path.isfile(pathOfAfile))
        os.remove("./test[EC:1.1.1.1].fasta")
        os.remove("./testHMMs/test[EC:1.1.1.1].aln")
        os.removedirs(autoMSAmafft.newDirectoryName)

if __name__ == "__main__":
    unittest.main()