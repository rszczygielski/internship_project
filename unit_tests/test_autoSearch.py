import sys
import unittest
import os
sys.path.insert(1, '/home/radeksz/bioinf/uni-gem-rec/ecnumModel/HMMaproach')
from unittest import TestCase
from autoSearch import HMMsearch

class AutoHmmTest(TestCase):

    def testAddFilesToList(self):
        with open("test.hmm", "w") as fileTest:
            pass
        with open("targetGenom.fasta", "w") as targetFile:
            targetFile.write(""">sp|P00968|CARB_ECOLI
MPKRTDIKSILILGAGPIVIGQACEFDYSGAQACKALREEGYRVILVNSNPATIMTDPEM
ADATYIEPIHWEVVRKIIEKERPDAVLPTMGGQTALNCALELERQGVLEEFGVTMIGATA
DAIDKAEDRRRFDVAMKKIGLETARSGIAHTMEEALAVAADVGFPCIIRPSFTMGGSGGG
IAYNREEFEEICARGLDLSPTKELLIDESLIGWKEYEMEVVRDKNDNCIIVCSIENFDAM
GIHTGDSITVAPAQTLTDKEYQIMRNASMAVLREIGVETGGSNVQFAVNPKNGRLIVIEM
NPRVSRSSALASKATGFPIAKVAAKLAVGYTLDELMNDITGGRTPASFEPSIDYVVTKIP
""") 
        cwd = os.getcwd()
        autoSearch = HMMsearch(cwd, "./targetGenom.fasta")
        pathOfAfile = os.path.join(cwd, "test.hmm")
        self.assertEqual(len(autoSearch.fileList), 1)
        self.assertEqual(autoSearch.fileList[0], pathOfAfile)
        os.remove("./test.hmm")
        os.remove("./targetGenom.fasta")
        os.removedirs(autoSearch.pathOfNewDir)
    
    def testCreateNewDirectory(self):
        with open("test.hmm", "w") as fileTest:
            pass
        with open("targetGenom.fasta", "w") as targetFile:
            targetFile.write(""">sp|P00968|CARB_ECOLI
MPKRTDIKSILILGAGPIVIGQACEFDYSGAQACKALREEGYRVILVNSNPATIMTDPEM
ADATYIEPIHWEVVRKIIEKERPDAVLPTMGGQTALNCALELERQGVLEEFGVTMIGATA
DAIDKAEDRRRFDVAMKKIGLETARSGIAHTMEEALAVAADVGFPCIIRPSFTMGGSGGG
IAYNREEFEEICARGLDLSPTKELLIDESLIGWKEYEMEVVRDKNDNCIIVCSIENFDAM
GIHTGDSITVAPAQTLTDKEYQIMRNASMAVLREIGVETGGSNVQFAVNPKNGRLIVIEM
NPRVSRSSALASKATGFPIAKVAAKLAVGYTLDELMNDITGGRTPASFEPSIDYVVTKIP
""") 
        cwd = os.getcwd()
        autoSearch = HMMsearch(cwd, "./targetGenom.fasta")
        self.assertTrue(os.path.isdir(autoSearch.pathOfNewDir))
        os.remove("./test.hmm")
        os.remove("./targetGenom.fasta")
        os.removedirs(autoSearch.pathOfNewDir)
    
    def testCreateHMMs(self):
        with open("test[EC:1.1.1.1].hmm", "w") as fileTest:
            fileTest.write("""HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y   
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   2.36423  4.10036  3.01026  2.63513  3.48914  2.60577  3.65657  2.80172  2.76527  2.51815  3.62117  3.25066  3.34819  3.14141  3.01880  2.75466  2.83338  2.43196  4.93243  3.74194
          2.68619  4.42223  2.77518  2.73124  3.46354  2.40513  3.72495  3.29355  2.67741  2.69355  4.24677  2.90347  2.73740  3.18147  2.89795  2.37887  2.77520  2.98519  4.58477  3.61504
          0.21067  3.15456  1.91530  0.74116  0.64733  0.00000        *
      1   2.90393  4.47562  3.98274  3.55268  3.17223  3.80827  4.30260  2.33867  3.31143  1.77849  1.55541  3.83947  4.25711  3.70203
""")
        with open("targetGenom.fasta", "w") as targetFile:
            targetFile.write(""">sp|P00968|CARB_ECOLI 
MPKRTDIKSILILGAGPIVIGQACEFDYSGAQACKALREEGYRVILVNSNPATIMTDPEM
ADATYIEPIHWEVVRKIIEKERPDAVLPTMGGQTALNCALELERQGVLEEFGVTMIGATA
DAIDKAEDRRRFDVAMKKIGLETARSGIAHTMEEALAVAADVGFPCIIRPSFTMGGSGGG
IAYNREEFEEICARGLDLSPTKELLIDESLIGWKEYEMEVVRDKNDNCIIVCSIENFDAM
GIHTGDSITVAPAQTLTDKEYQIMRNASMAVLREIGVETGGSNVQFAVNPKNGRLIVIEM
NPRVSRSSALASKATGFPIAKVAAKLAVGYTLDELMNDITGGRTPASFEPSIDYVVTKIP
""") 
        cwd = os.getcwd()
        autoSearch = HMMsearch(cwd, "./targetGenom.fasta")
        self.assertTrue(os.path.isdir(autoSearch.pathOfNewDir))
        autoSearch.searchTargetGenome(os.cpu_count())
        pathOfAfile = os.path.join(cwd, "./searchDir/test[EC:1.1.1.1].out")
        self.assertTrue(os.path.isfile(pathOfAfile))
        os.remove("./test[EC:1.1.1.1].hmm")
        os.remove("./searchDir/test[EC:1.1.1.1].out")
        os.remove("./targetGenom.fasta")
        os.removedirs(autoSearch.pathOfNewDir)

if __name__ == "__main__":
    unittest.main()