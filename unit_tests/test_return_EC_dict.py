import unittest
import sys
import os
sys.path.insert(1, '/home/radeksz/bioinf/uni-gem-rec/ecnumModel')
from unittest import TestCase
import return_EC_dict

class SaveEcNumsTest(TestCase):
    
    def testDeleteBreackets(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            pass
        self.assertTrue(os.path.isfile(testPath))
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        testList = ["PMK1_01696(adhB)", "PMK1_04035(adh)", "VK055_0781 VK055_3271(fucO2)", "VK055_0781", "VK055_3271"]
        testList = return_EC.deleteBrackets(testList)
        for elem in testList:
            self.assertNotIn("(", elem)
            self.assertNotIn(")", elem)
        os.remove(testPath)
        
    
    def testReturningListIfManyECnumsInLine(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            pass
        self.assertTrue(os.path.isfile(testPath))
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        os.remove(testPath)
        testEC = "[EC:1.1.1.4 1.1.1.- 1.1.1.303 1.2.3.403]"
        testListECnums = return_EC.returnListOfECnumsIfmanyECnumsInLine(testEC)
        correctList = ["[EC:1.1.1.4]", "[EC:1.1.1.303]", "[EC:1.2.3.403]"]
        self.assertListEqual(testListECnums, correctList)

    def testReturningListIfManyECnumsInLineWithGap(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            pass
        self.assertTrue(os.path.isfile(testPath))
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        os.remove(testPath)
        testEC = "[EC:1.1.1.4 1.1.1.- 1.1.1.303 1.2.3.-]"
        testListECnums = return_EC.returnListOfECnumsIfmanyECnumsInLine(testEC)
        correctList = ["[EC:1.1.1.4]", "[EC:1.1.1.303]"]
        self.assertListEqual(testListECnums, correctList)

    def testReturnECnumIfECnumLine(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            pass
        self.assertTrue(os.path.isfile(testPath))
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        testECline = "DEFINITION  alcohol dehydrogenase [EC:1.1.1.1]"
        testEC = return_EC.returnECnumIfECnumLine(testECline)
        correctEC = "[EC:1.1.1.1]"
        self.assertEqual(testEC, correctEC)
        os.remove(testPath)
    
    def testReturnECnumIfECnumLineWithGap(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            pass
        self.assertTrue(os.path.isfile(testPath))
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        testECline = "DEFINITION  alcohol dehydrogenase [EC:1.1.1.-]"
        testEC = return_EC.returnECnumIfECnumLine(testECline)
        correctEC = None
        self.assertEqual(testEC, correctEC)
        os.remove(testPath)
        
    def testSaveECnumToFileWithoutID(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            testFile.write("ENTRY       K00001                      KO\n")
            testFile.write("DEFINITION  alcohol dehydrogenase [EC:1.1.1.1]\n")
            testFile.write("///")
        self.assertTrue(os.path.isfile(testPath))
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        return_EC.returnEcIDsDict()
        return_EC.saveECnumToFile()
        with open(return_EC.path) as saveFile:
            lineWithEC = saveFile.readlines()[0]
        self.assertAlmostEqual(lineWithEC.strip(), "K00001-[EC:1.1.1.1]")
        os.remove(return_EC.path)
        os.remove(testPath)
    
    def testSaveECnumToFileWithID(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            testFile.write("ENTRY       K00001                      KO\n")
            testFile.write("DEFINITION  alcohol dehydrogenase [EC:1.1.1.1]\n")
            testFile.write("GENES       DME: Dmel_CG3481(Adh)\n")
            testFile.write("\t\t\tDMO: Dmoj_GI17643(Dmoj_Adh2) Dmoj_GI17644(Dmoj_Adh1)\n")
            testFile.write("\t\t\tDYA: Dyak_GE19037(Dyak_Adh)\n")
            testFile.write("///")
        self.assertTrue(os.path.isfile(testPath))
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        return_EC.returnEcIDsDict()
        return_EC.saveECnumToFile(withIDs=True)
        with open(return_EC.path) as saveFile:
            listOfLines = saveFile.readlines()
            lineWithEC = listOfLines[0]
            secondLine = listOfLines[1]
            thirdLine = listOfLines[2]
            fourthLine = listOfLines[3]
            fifthLine = listOfLines[4]
        self.assertAlmostEqual(lineWithEC.strip(), "K00001-[EC:1.1.1.1]:")
        self.assertAlmostEqual(secondLine.strip(), "Dmel_CG3481")
        self.assertAlmostEqual(thirdLine.strip(), "Dmoj_GI17643")
        self.assertAlmostEqual(fourthLine.strip(), "Dmoj_GI17644")
        self.assertAlmostEqual(fifthLine.strip(), "Dyak_GE19037")
        os.remove(return_EC.path)
        os.remove(testPath)

    def testAddKeyValueToDictAndResetVariablesEcNotList(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            pass
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        return_EC.ecNum = "[EC:1.1.1.1]"
        return_EC.koNum = "K00001"
        return_EC.peptideIDslist = ["Dmel_CG3481", "Dmoj_GI17643", "Dmoj_GI17644"]
        return_EC.addKeyValueToDictAndResetVariables()
        self.assertEqual(return_EC.ecNumIDsDict["K00001-[EC:1.1.1.1]"], ["Dmel_CG3481", "Dmoj_GI17643", "Dmoj_GI17644"])
        os.remove(testPath)
    
    def testAddKeyValueToDictAndResetVariablesEcIsAList(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            pass
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        return_EC.ecNum = ["[EC:1.1.1.1]", "[EC:1.1.1.2]"]
        return_EC.koNum = "K00001"
        return_EC.peptideIDslist = ["Dmel_CG3481", "Dmoj_GI17643", "Dmoj_GI17644"]
        return_EC.idBool = True
        return_EC.addKeyValueToDictAndResetVariables()
        self.assertEqual(return_EC.ecNumIDsDict["K00001-[EC:1.1.1.1]"], ["Dmel_CG3481", "Dmoj_GI17643", "Dmoj_GI17644"])
        self.assertEqual(return_EC.ecNumIDsDict["K00001-[EC:1.1.1.2]"], ["Dmel_CG3481", "Dmoj_GI17643", "Dmoj_GI17644"])
        self.assertFalse(return_EC.idBool)
        self.assertEqual(return_EC.ecNum, "")
        self.assertEqual(return_EC.koNum, "")
        self.assertEqual(return_EC.peptideIDslist, [])
        os.remove(testPath)
    
    def testChangeIdBoolIntoTrue(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            pass
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        self.assertFalse(return_EC.idBool)
        return_EC.changeIdBoolIntoTrue("GENES       DME: Dmel_CG3481(Adh)\n".split())
        self.assertTrue(return_EC.idBool)
        self.assertEqual(return_EC.peptideIDslist, ["Dmel_CG3481"])
        os.remove(testPath)

    def testReturnEcNumIDsDict(self):
        cwd = os.getcwd()
        testPath = os.path.join(cwd, 'testFastaFile.fasta')
        self.assertFalse(os.path.isfile(testPath))
        with open(testPath, "w") as testFile:
            testFile.write("ENTRY       K00001                      KO\n")
            testFile.write("DEFINITION  alcohol dehydrogenase [EC:1.1.1.1]\n")
            testFile.write("GENES       DME: Dmel_CG3481(Adh)\n")
            testFile.write("\t\t\tDMO: Dmoj_GI17643(Dmoj_Adh2) Dmoj_GI17644(Dmoj_Adh1)\n")
            testFile.write("\t\t\tDYA: Dyak_GE19037(Dyak_Adh)\n")
            testFile.write("///")
        self.assertTrue(os.path.isfile(testPath))
        return_EC = return_EC_dict.ReturnEcIDs(testPath)
        return_EC.returnEcIDsDict()
        print(return_EC.ecNumIDsDict)
        self.assertEqual(return_EC.ecNumIDsDict["K00001-[EC:1.1.1.1]"], ["Dmel_CG3481","Dmoj_GI17643", "Dmoj_GI17644", "Dyak_GE19037"])
        os.remove(testPath)

if __name__ == "__main__":
    unittest.main()
    # Dmel_CG3481 
	# Dmoj_GI17643 
	# Dmoj_GI17644 
	# Dyak_GE19037 