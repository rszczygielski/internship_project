import os
import sys
import unittest
from unittest import TestCase
sys.path.insert(1, '/home/radeksz/bioinf/uni-gem-rec/classificationDir')
from classification import Classification

class ClassificationTest(TestCase):

    @classmethod
    def setUpClass(cls):
        with open("uniprotTestFile.tab", "w") as uniprotFile:
            uniprotFile.write("""Entry	EC number
P0AE70	3.1.27.-
P0ACF8
P0ABG4	2.4.1.129
P23837	2.7.13.3; 3.1.3.-
P23874	2.7.11.1
P76129	3.1.4.52
P0A9V8	1.1.1.373; 1.1.1.61
""")
        with open("searchTestFile.txt", "w") as searchFile:
            searchFile.write("""Query ID	Predicted EC number
P0ABG4  2.4.1.129
P23837  2.7.13.3
P23874  2.7.11.1; 1.1.1.61
P76129  3.1.4.52
P0A9V8  1.1.1.373; 1.1.1.61
P0ACF0  1.1.1.61""")
        
        with open("deepecTest.txt", "w") as deepResult:
            deepResult.write("""Query ID	Predicted EC number
sp|P04982|RBSD_ECOLI	EC:5.4.99.62
sp|P04994|EX7L_ECOLI	EC:3.1.11.6
sp|P0A7Z4|RPOA_ECOLI	EC:2.7.7.6
sp|P0A962|ASPG1_ECOLI	EC:3.5.1.1
sp|P0A9Q7|ADHE_ECOLI	EC:1.2.1.10
sp|P0A9Q7|ADHE_ECOLI	EC:1.1.1.1
sp|P0ABF8|PGSA_ECOLI	EC:2.7.8.5""")
        cwd = os.getcwd()
        cls.uniProtTestFilePath = os.path.join(cwd, "uniprotTestFile.tab")
        cls.searchFilePath = os.path.join(cwd, "searchTestFile.txt")
        cls.deepEcFilePath = os.path.join(cwd, "deepecTest.txt")

    @classmethod
    def tearDownClass(cls):
        os.remove(__class__.uniProtTestFilePath)
        os.remove(__class__.searchFilePath)
        os.remove(__class__.deepEcFilePath) 


    def testCreateUniprotDictWithoutCutBool(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=False, deepBool=False)
        expectedDict = {'P0ABG4': ['2.4.1.129'], 'P23837': ['2.7.13.3'], 
        'P23874': ['2.7.11.1'], 'P76129': ['3.1.4.52'], 'P0A9V8': ['1.1.1.373', '1.1.1.61']}
        self.assertEqual(expectedDict, classification.uniProtDict)
        self.assertEqual(len(classification.uniProtDict), 5)
        
    
    def testCreateUniprotDictWithCutBool(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=True, deepBool=False)
        expectedDict = {'P0ABG4': ['2.4.1'], 'P23837': ['2.7.13'], 'P23874': ['2.7.11'],
        'P76129': ['3.1.4'], 'P0A9V8': ['1.1.1', '1.1.1']}
        self.assertEqual(expectedDict, classification.uniProtDict)
        self.assertEqual(len(classification.uniProtDict), 5) 
    
    def testCreateSearchDictWithoutCutBool(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=False, deepBool=False)
        expectedDict = {'P0ABG4': ['2.4.1.129'], 'P23837': ['2.7.13.3'], 
        'P23874': ['2.7.11.1', "1.1.1.61"], 'P76129': ['3.1.4.52'], 'P0A9V8': ['1.1.1.373', '1.1.1.61'], "P0ACF0": ["1.1.1.61"]}
        self.assertEqual(expectedDict, classification.searchDict)
        self.assertEqual(len(classification.searchDict), 6)

    def testCreateSearchDictWithCutBool(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=True, deepBool=False)
        expectedDict = {'P0ABG4': ['2.4.1'], 'P23837': ['2.7.13'], 
        'P23874': ['2.7.11', "1.1.1"], 'P76129': ['3.1.4'], 'P0A9V8': ['1.1.1', '1.1.1'], "P0ACF0": ["1.1.1"]}
        self.assertEqual(expectedDict, classification.searchDict)
        self.assertEqual(len(classification.searchDict), 6) 
    
    def testCreateDeepDict(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.deepEcFilePath, cutBool=False, deepBool=True)
        expectedDict = {'P04982': ['5.4.99.62'], 'P04994': ['3.1.11.6'], 'P0A7Z4': ['2.7.7.6'],
         'P0A962': ['3.5.1.1'], 'P0A9Q7': ['1.2.1.10', '1.1.1.1'], 'P0ABF8': ['2.7.8.5']}
        self.assertEqual(expectedDict, classification.searchDict)
        self.assertEqual(len(classification.searchDict), 6)
    
    def testAddTruePositivesAndFalsePositives(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=False, deepBool=False)
        self.assertEqual(classification.tp, 6)
        self.assertEqual(classification.fp, 1)
        classification.uniProtDict["test1"] = ["uniProtEC1"]
        classification.uniProtDict["test2"] = ["uniProtEC2"]
        classification.searchDict["test1"] = ["wrong1"]
        classification.searchDict["test2"] = ["wrong2", "uniProtEC2"]
        classification.fp = 0
        classification.tp = 0
        classification.addTruePositivesAndFalsePositives()
        self.assertEqual(classification.fp, 3)
        self.assertEqual(classification.tp, 7)


    
    def testAddFalseNegatives(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=False, deepBool=False)
        classification.uniProtDict["test1"] = ["uniProtEC1"]
        classification.uniProtDict["test2"] = ["uniProtEC2"]
        classification.searchDict["test1"] = ["searchEC1"]
        classification.searchDict["test2"] = ["searchEC2"]
        self.assertEqual(classification.fn, 0)
        classification.addFalseNegatives()
        self.assertEqual(classification.fn, 2)
        self.assertNotEqual(classification.searchDict, classification.uniProtDict)

    def testCutECnum(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=False, deepBool=False)
        testECNum = "1.1.1.61"
        ecNumAfterCut = classification.cutEcNum(testECNum)
        self.assertEqual(ecNumAfterCut, "1.1.1")
    
    def testReturnPrecision(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=False, deepBool=False)
        zeroPrecision = classification.returnPrecision(tp=0, fp=100)
        self.assertEqual(zeroPrecision, 0)
        fiftyPrecision = classification.returnPrecision(tp=50, fp=50)
        self.assertEqual(fiftyPrecision, 0.5)

    def testReturnRecall(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=False, deepBool=False)
        zeroRecall = classification.returnRecall(tp=0, fn=100)
        self.assertEqual(zeroRecall, 0)
        fiftyRecall = classification.returnRecall(tp=50, fn=50)
        self.assertEqual(fiftyRecall, 0.5)
    
    def testJaccardIndex(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=False, deepBool=False)
        del classification.searchDict["P0ACF0"]
        classification.searchDict["P23874"].remove("1.1.1.61")
        self.assertDictEqual(classification.searchDict, classification.uniProtDict)
        jaccardIndex = float(classification.jaccardIndex().split()[-1])
        self.assertEqual(jaccardIndex, 1.0)
        for i in range(len(set(classification.uniProtEcNumList + classification.searchEcNumList))):
            classification.searchEcNumList.append(i)        
        jaccardIndex = float(classification.jaccardIndex().split()[-1])
        self.assertEqual(jaccardIndex, 0.5)
    
    def testAddTrueNegatives(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=False, deepBool=False)
        self.assertEqual(classification.tn, 1)
        classification.emptyUniProtPeptideIDs.append("test")
        classification.tn = 0
        classification.addTrueNegatives()
        self.assertEqual(classification.tn, 2)
    
    def testReturnSpecificity(self):
        classification =  Classification(self.__class__.uniProtTestFilePath, self.__class__.searchFilePath, cutBool=False, deepBool=False)
        self.assertEqual(classification.tn, 1)
        self.assertEqual(classification.fp, 1)
        print(classification.tn)
        self.assertEqual(classification.returnSpecificity(classification.tn, classification.fp), 0.5)
        self.assertEqual(classification.returnSpecificity(1, 1), 0.5)
        self.assertEqual(classification.returnSpecificity(tn=100, fp=50), 0.67)
        classification.uniProtDict["test1"] = ["uniProtEC1"]
        classification.uniProtDict["test2"] = ["uniProtEC2"]
        classification.searchDict["test1"] = ["wrong1"]
        classification.searchDict["test2"] = ["wrong2"]
        classification.fp = 0
        classification.addTruePositivesAndFalsePositives()
        self.assertEqual(classification.returnSpecificity(classification.tn, classification.fp), 0.25)
        
        





if __name__ == "__main__":
    unittest.main()