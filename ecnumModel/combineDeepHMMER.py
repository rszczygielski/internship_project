import os
import argparse
import pandas as pd

class CombineDeepHmmer():
    def __init__(self, deepResult, hmmerResult, fileName):
        self.deepResult = deepResult
        self.hmmerResult = hmmerResult
        self.fileName = fileName

    def returnHmmerDict(self, hmmerResult):
        dataSereisDict = dict(pd.read_table(hmmerResult, usecols=["Query ID", "Predicted EC number"], index_col="Query ID").squeeze('columns'))
        for key, value in dataSereisDict.items():
            splitedValue = value.split()
            dataSereisDict[key] = splitedValue
        return dataSereisDict


    def returnDeepDict(self, deepResult):
        deepDict = {}
        with open(deepResult) as deepFile:
            for line in deepFile.readlines()[1:]:
                # loop through HMMER search file lines  
                splitedLine = line.strip().split()
                # print(splitedLine)
                splitedID = splitedLine[0].split("|")
                peptideID = splitedID[1]
                searchEcNum = splitedLine[1][3:]
                if "n" in searchEcNum:
                        continue
                if peptideID not in deepDict:
                    deepDict[peptideID] = [searchEcNum]
                else:
                    deepDict[peptideID].append(searchEcNum)
        return deepDict
    
    def returnCombinedDict(self, deepDict, hmmerDict):
        combinedDict = deepDict
        for peptideID in hmmerDict:
            # print(hmmerDict[peptideID])
            if peptideID in combinedDict:
                for ecNumber in hmmerDict[peptideID]:
                    # print(type(hmmerDict[peptideID]))
                    if ecNumber not in combinedDict[peptideID]:
                        combinedDict[peptideID].append(ecNumber)
                    else:
                        continue
            else:
                combinedDict[peptideID] = hmmerDict[peptideID]
        return combinedDict
    
    def saveEcKoToFile(self):
        deepDict = self.returnDeepDict(self.deepResult)
        hmmerDict = self.returnHmmerDict(self.hmmerResult)
        combinedDict = self.returnCombinedDict(deepDict, hmmerDict)
        with open(self.fileName, "w") as saveFile:
            saveFile.write("Query ID\tPredicted EC number\n")
            for peptideID in combinedDict:
                listToStringOfECs = " ".join(combinedDict[peptideID])
                # listToStringOfECs = list(combinedDict[peptideID])
                saveFile.write(f"{peptideID}\t{listToStringOfECs}\n")

if __name__ == "__main__":
    combineDeepHmmer = CombineDeepHmmer("/home/radeksz/bioinf/essential_files/DeepAproach/E.coliDeep/DeepEC_Result.txt", "/home/radeksz/bioinf/uni-gem-rec/ECnumsAfterTshBestHit*-10^150.txt", "CombinedECs")
    combineDeepHmmer.saveEcKoToFile()