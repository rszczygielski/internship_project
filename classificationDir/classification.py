import argparse

class Classification():
    def __init__(self, uniprotFile, searchFile, cutBool, deepBool):
        self.uniprotFile = uniprotFile
        self.searchFile = searchFile
        self.cutBool = cutBool
        self.tp = 0
        self.fp = 0
        self.fn = 0
        self.tn = 0
        self.searchEcNumList = []
        self.uniProtEcNumList = []
        self.emptyUniProtPeptideIDs = []
        self.searchDict = {}
        self.uniProtDict = {}
        self.createUniprotDict(cutBool)
        if deepBool:
            self.createDeepDict(cutBool)
        else:
            self.createSearchDict(cutBool)
        self.addTruePositivesAndFalsePositives()
        self.addFalseNegatives()
        print(len(self.uniProtDict))
        # self.addTrueNegatives()
        
    def createUniprotDict(self, cutBool):
        """Method creates dictionary with peptide ID as a key and EC numbers list as a value based on uniProt file
        """
        with open(self.uniprotFile) as uniprotFile:
            for line in uniprotFile.readlines()[1:]:
                splitedLine = line.split()
                peptideID = splitedLine[0]
                if len(splitedLine) == 1:
                    self.emptyUniProtPeptideIDs.append(peptideID)
                    continue
                ecNumList = []
                ecNums = splitedLine[1:]
                for ecNum in ecNums:
                    if "-" in ecNum or "n" in ecNum:
                        continue
                    if ";" in ecNum:
                        ecNum = ecNum[:-1]
                    if cutBool:
                        ecNum = self.cutEcNum(ecNum)
                    ecNumList.append(ecNum)
                    if ecNum not in self.uniProtEcNumList:
                        self.uniProtEcNumList.append(ecNum)
                if len(ecNumList) != 0:
                    self.uniProtDict[peptideID] = ecNumList

    def createDeepDict(self, cutBool):
        """Method creates a dictionary where key is a peptide ID from search and value is a list of estimated EC nums using deep approach
        """
        with open(self.searchFile) as searchFile:
            for line in searchFile.readlines()[1:]:
                # loop through HMMER search file lines  
                splitedLine = line.strip().split()
                # print(splitedLine)
                splitedID = splitedLine[0].split("|")
                peptideID = splitedID[1]
                searchEcNum = splitedLine[1][3:]
                if "-" in searchEcNum or "n" in searchEcNum:
                        continue
                if cutBool:
                    searchEcNum = self.cutEcNum(searchEcNum)
                if searchEcNum not in self.searchEcNumList:
                        self.searchEcNumList.append(searchEcNum)
                if peptideID not in self.searchDict:
                    self.searchDict[peptideID] = [searchEcNum]
                else:
                    self.searchDict[peptideID].append(searchEcNum)

    def createSearchDict(self, cutBool):
        with open(self.searchFile) as searchFile:
            for line in searchFile.readlines()[1:]:
                splitedLine = line.split()
                peptideID = splitedLine[0]
                ecNumList = []
                ecNums = splitedLine[1:]
                for ecNum in ecNums:
                    if "-" in ecNum or "n" in ecNum:
                        continue
                    if ";" in ecNum:
                        ecNum = ecNum[:-1]
                    if cutBool:
                        ecNum = self.cutEcNum(ecNum)
                    if ecNum not in self.searchEcNumList:
                        self.searchEcNumList.append(ecNum)
                    ecNumList.append(ecNum)
                if len(ecNumList) != 0:
                    self.searchDict[peptideID] = ecNumList
                
    def addTruePositivesAndFalsePositives(self):
        """Method adds a True positives and False positives if it find one
        """
        for searchPeptideID in self.searchDict:
            if searchPeptideID in self.uniProtDict:
                for ecNumber in self.searchDict[searchPeptideID]:
                    if ecNumber in self.uniProtDict[searchPeptideID]:
                        self.tp += 1
                    if ecNumber not in self.uniProtDict[searchPeptideID]:
                        self.fp += 1
    
    def addFalseNegatives(self):
        """Method adds a False negatives if it find one
        """
        for uniProtPeptideID in self.uniProtDict:
            if uniProtPeptideID in self.searchDict:
                for ecNumber in self.uniProtDict[uniProtPeptideID]:
                    if ecNumber not in self.searchDict[uniProtPeptideID]:
                        self.fn += 1
    
    def addTrueNegatives(self):
        for uniProtPeptideID in self.emptyUniProtPeptideIDs:
            if uniProtPeptideID not in self.searchDict:
                self.tn += 1
    
    def cutEcNum(self, ecNum):
        """Method cuts a last digit of a EC number and returns it

        Args:
            ecNum (str): EC number

        Returns:
            str: EC without last digit
        """
        splitedEcNum = ecNum.split(".")
        listWithoutLastNum = splitedEcNum[:-1]
        newEcNum = ".".join(listWithoutLastNum)
        return newEcNum

    def returnPrecision(self, tp, fp):
        """Method returns a precision score

        Args:
            tp (int): true positive
            fp (int): false positive

        Returns:
            int: Precision score
        """
        if tp == 0:
            return 0
        else:
            return round(tp/(tp+fp), 2)
    
    def returnRecall(sefl, tp, fn):
        """Method returns a recall score

        Args:
            tp (int): true positive
            fn (int): false negative

        Returns:
            int: Recall score
        """
        if tp == 0:
            return 0
        else:
            return round(tp/(tp+fn), 2)
    
    def returnSpecificity(self, tn, fp):
        """Method returns a specificity score

        Args:
            tn (int): true negative
            fn (int): false negative

        Returns:
            int: Specificity score
        """
        if tn == 0:
            return 0
        else:
            return round(tn/(tn+fp), 2)
    
    def jaccardIndex(self):
        """Method counts a Jaccard Index and prints it
        """
        intersection = 0
        union = len(set(self.uniProtEcNumList + self.searchEcNumList))
        for ecNumber in self.searchEcNumList:
            if ecNumber in self.uniProtEcNumList:
                intersection += 1
        jaccardIndexValue = intersection/union
        return f"Intersection: {intersection}\nUnion: {union}\nJaccard Index: {round(jaccardIndexValue, 3)}"
        
    def saveOverallReacallAndPrecision(self, outputName):
        with open(outputName, "w") as saveFile:
            overallPrecision = self.returnPrecision(self.tp, self.fp)
            overallRecall = self.returnRecall(self.tp, self.fn)
            # overallSpecificity = self.returnSpecificity(self.tn, self.fp)
            hashLine = "#"*42
            saveFile.write(f"{hashLine}\n")
            saveFile.write(f"OVERALL PRECISION SCORE: {overallPrecision}\n")
            saveFile.write(f"{hashLine}\n")
            saveFile.write(f"OVERALL RECALL SCORE: {overallRecall}\n")
            saveFile.write(f"{hashLine}\n")
            # saveFile.write(f"OVERALL Specificity SCORE: {overallSpecificity}\n")
            # saveFile.write(f"{hashLine}\n")
            saveFile.write(f"ALL TRUE POSITIVES {self.tp}\n")
            # saveFile.write(f"ALL TRUE NEGATIVES {self.tn}\n")
            saveFile.write(f"ALL FALSE POSITIVES {self.fp}\n")
            saveFile.write(f"ALL FALSE NEGATIVES {self.fn}\n")
            saveFile.write(f"{hashLine}\n")
            saveFile.write(f"{self.jaccardIndex()}\n")
            saveFile.write(f"{hashLine}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Program runs a calassification for results form previous runs")
    parser.add_argument("reference", metavar="reference", help="Path to referance for classification")
    parser.add_argument("results", metavar="results",help="Path to file with results")
    parser.add_argument("-o", metavar="output" , dest="output",default="Score.txt" ,help="change output name -> default Score.txt")
    parser.add_argument("-c", dest="cut", action="store_true", help="cut the last digit from EC number - default False")
    parser.add_argument("-d", dest="deep", action="store_true", help="change classification for deep approach, analise results from deepec")
    args = parser.parse_args()
    reference = args.reference
    results = args.results
    output = args.output
    cut = args.cut
    deep_bool = args.deep
    
    classification = Classification(reference, results, cut, deep_bool)
    classification.saveOverallReacallAndPrecision(output)

# reference = /home/radeksz/bioinf/uni-gem-rec/classificationDir/uniprot-organism_escherichia+coli+k12-filtered-reviewed_yes.tab