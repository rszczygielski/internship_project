import re
import pandas as pd
import argparse
import os
from logger import Logger

class MNXrefID():
    def __init__(self, mnxID, mnxEquation, reference, isBalanced, isTransport, peptideID, ecNumber):
        self.mnxID = mnxID
        self.mnxEquation = mnxEquation
        self.reference = reference
        self.isBalanced = isBalanced
        self.isTransport = isTransport
        self.peptideID = peptideID
        self.ecNumber = ecNumber

class Reactoins():
    def __init__(self, mnxRefFile, ecNumResult):
        self.mnxRefFile = mnxRefFile
        self.ecNumResult = ecNumResult
        self.temporaryFile = "dataFrameNotNormalize.txt"
        self.mnxRefIDstrucList = []
        Logger.settings(show_date=False, show_file_name=False)
        self.addAllPossibleMNXrefStructsToList()
    
    def removeFileDecorator(method):
        def wrapper(self, *args):
            result = method(self, *args)
            self.deleteFile()
            return result
        return wrapper

    def readMnxRefFile(self):
        mnxDataFrame = pd.read_table(self.mnxRefFile, header=351, usecols=["ID", 
        "mnx_equation", "reference", "classifs", "is_balanced",	"is_transport"])
        mnxDataFrameOnlyBalanced = mnxDataFrame.loc[mnxDataFrame["is_balanced"] == "B"]
        return mnxDataFrameOnlyBalanced[mnxDataFrameOnlyBalanced["classifs"].notna()]

    def readEcNumResultFileAndReturnDict(self):
        ecNumResultDataFrame = pd.read_table(self.ecNumResult, usecols=["Query ID", "Predicted EC number"],
        index_col="Query ID").squeeze("columns")
        ecNumResultDict = dict(ecNumResultDataFrame)
        for peptideID, ecNumber in ecNumResultDict.items():
            ecNumberSplited = ecNumber.split()
            ecNumResultDict[peptideID] = ecNumberSplited
        return ecNumResultDict
    
    def replaceEqualSighFromEquation(self, equation):
        equation = re.sub(r"@(\w+)", r"[\1]", equation)
        return equation.replace("=", "<=>")
    
    def addOneMNXidStructToList(self, ecNumber, peptideID, mnxDataFrame):
        rowWithMetabolite = mnxDataFrame.loc[mnxDataFrame["ID"] == "MNXM1000000"]
        numpyArrayRowsListed = rowWithMetabolite.values.astype(str)
        for row in numpyArrayRowsListed:
            mnxID = row[0]
            mnxEquation = self.replaceEqualSighFromEquation(row[1])
            reference = row[2]
            isBalanced = row[3]
            isTransport = row [4]
            self.mnxRefIDstrucList.append(MNXrefID(mnxID, mnxEquation, reference, isBalanced, isTransport, peptideID, ecNumber))
            Logger.INFO(f"MNXref struct with ID: {mnxID} added to the list")

    def addAllPossibleMNXrefStructsToList(self):
        ecNumResultDict = self.readEcNumResultFileAndReturnDict()
        mnxDataFrame = self.readMnxRefFile()
        for peptideID, ecNumberList in ecNumResultDict.items():
            for ecNumber in ecNumberList:
                self.addOneMNXidStructToList(ecNumber,peptideID, mnxDataFrame)

    @removeFileDecorator
    def reorderMNXrefIDs(self, saveFile):
        mnxDataFrame = pd.read_table(self.temporaryFile)
        uniqueListOfMNXrefIDs = list(mnxDataFrame["ID"].unique())
        for mnxID in uniqueListOfMNXrefIDs:
            theSameMnxIDs = mnxDataFrame.loc[mnxDataFrame["ID"] == mnxID]
            listOfpeptideIDs = list(theSameMnxIDs["peptideID"].unique())
            listOfEcNumbers = list(theSameMnxIDs["EC"].unique())
            indexesOftheSameMnxID = theSameMnxIDs.index
            theSameMnxIDs = theSameMnxIDs.assign(peptideID = " ".join(listOfpeptideIDs), EC = " ".join(listOfEcNumbers))
            oneUniqueMnxID = theSameMnxIDs.drop(indexesOftheSameMnxID[1:])
            mnxDataFrame = mnxDataFrame.drop(indexesOftheSameMnxID)
            mnxDataFrame = pd.concat([mnxDataFrame, oneUniqueMnxID])
        mnxDataFrame.to_csv(saveFile, sep='\t', index=False)

    def saveResultsToFile(self):
        with open(self.temporaryFile, "w") as resultFile:
            resultFile.write("ID\tmnx_equation\tEC\tpeptideID\treference\tis_balanced\tis_transport\n")
            for mnxRefID in self.mnxRefIDstrucList:
                resultFile.write(f"{mnxRefID.mnxID}\t{mnxRefID.mnxEquation}\t{mnxRefID.ecNumber}\t{mnxRefID.peptideID}\t{mnxRefID.reference}\t{mnxRefID.isBalanced}\t{mnxRefID.isTransport}\n")

    def deleteFile(self):
        os.remove(self.temporaryFile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map results with the MNXref data base referance")
    parser.add_argument("mnxRefFile", metavar="mnxRefFile", help="Path to MNXref data base file")
    parser.add_argument("ecNumResult", metavar="ecNumResult", help="Path to EC number result file")
    parser.add_argument("-o", metavar="output", dest="output", default="All_MNX_IDs_Combined.txt",
    help="Name of the output file. Default -> All_MNX_IDs_Combined.txt" )
    args = parser.parse_args()
    mnxRefFile = args.mnxRefFile
    ecNumResult = args.ecNumResult
    output = args.output
    reactoins = Reactoins(mnxRefFile, ecNumResult)
    reactoins.saveResultsToFile()
    reactoins.reorderMNXrefIDs(output)
    # reactoins = Reactoins("/home/radeksz/Downloads/reac_prop.tsv", 
    # "/home/radeksz/bioinf/uni-gem-rec/CombinedECs")
    # print(reactoins.readMnxRefFile()["classifs"])
    # print(reactoins.readEcNumResultFileAndReturnDict())
    # print(reactoins.addMNXidStructToList("1.2.1.3"))
    # reactoins.saveResultsToFile("PeptideWithMNXlines.txt")
    # # print(reactoins.returnNewDataFrame())
    # reactoins.reorderMNXrefIDs("finalOutput.txt")
    