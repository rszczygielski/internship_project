import pandas as pd



class SetOfECs():
    def __init__(self, deepResult, hmmResult):
        self.deepResult = deepResult
        self.hmmResult = hmmResult


    def returnSetOfDeepRsult(self):
        dataStructure = pd.read_table(self.deepResult, usecols=["Predicted EC number"])
        dataStructureWithoutNotFullECs = dataStructure[dataStructure["Predicted EC number"].str.contains("n") == False]
        return set(dataStructureWithoutNotFullECs["Predicted EC number"].str.replace("EC:", "").unique())
    
    def returnSetOfHmmResult(self):
        dataStructure = pd.read_table(self.hmmResult, usecols=["Predicted EC number"])
        listOfListsECs = dataStructure["Predicted EC number"].str.split().tolist()
        # listOfECs = [row for row in listOfListsECs if isinstance(row, list) row]
        hmmListOfECs = []
        for row in listOfListsECs:
            hmmListOfECs.extend(row)
        return set(hmmListOfECs)

    def countSetDifference(self):
        hmmSetOfECs = self.returnSetOfHmmResult()
        deepSetOfECs = self.returnSetOfDeepRsult()
        intersection = set([hmmEC for hmmEC in hmmSetOfECs if hmmEC in deepSetOfECs])
        lenDeepSetDifference = len(deepSetOfECs.difference(intersection))
        lenHmmSetDifference = len(hmmSetOfECs.difference(intersection))
        hmmSetOfECs.update(deepSetOfECs)
        union = len(hmmSetOfECs)
        deepResult = round(lenDeepSetDifference / union, 2)
        hmmResult = round(lenHmmSetDifference / union, 2)
        return deepResult, hmmResult, union, len(intersection), lenDeepSetDifference, lenHmmSetDifference

    def saveResults(self, output):
        with open(output, "w") as resultFile:
            deepResult, hmmResult, union, lenIntersection,  lenDeepSetDifference, lenHmmSetDifference = self.countSetDifference()
            hashLine = "#" * 45
            resultFile.write(f"{hashLine}\n")
            resultFile.write(f"DEEP LEARING APPROACH SET DIFFERENCE: {deepResult}\n")
            resultFile.write(f"{hashLine}\n")
            resultFile.write(f"HMM APPROACH SET DIFFERENCE: {hmmResult}\n")
            resultFile.write(f"{hashLine}\n")
            resultFile.write(f"Union lenght: {union}\n")
            resultFile.write(f"Intersection lenght: {lenIntersection}\n")
            resultFile.write(f"Deep Set Difference lenght: {lenDeepSetDifference}\n")
            resultFile.write(f"HMM Set Difference lenght: {lenHmmSetDifference}\n")
            resultFile.write(f"{hashLine}\n")



        




if __name__ == "__main__":
    setOfECs = SetOfECs("/home/radeksz/bioinf/essential_files/DeepAproach/E.coliDeep/DeepEC_Result.txt", 
    "/home/radeksz/bioinf/uni-gem-rec/ECnumsAfterTshBestHit*-10^150.txt")
    setOfECs.saveResults("SetDifference.txt")