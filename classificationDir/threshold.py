import os
import argparse

class Threshold():
    def __init__(self, dirPath, fileName, threshold, oneHitBool):
        self.dirPath = dirPath
        self.listOfFilePaths = []
        self.ecIdAfterThreshDir = {}
        self.peptideIDecListDict = {}
        self.fileName = fileName
        self.threshold = threshold
        self.addFilesPathToList(".out",self.dirPath)
        self.addToDict(oneHitBool)
        
    def addFilesPathToList(self, ext, direcotryPath):
        """Method adds all fasta format files from directory to a list
        """
        for fileInDir in os.listdir(direcotryPath):
            if fileInDir.endswith(ext):
                filePath = os.path.join(direcotryPath, fileInDir)
                self.listOfFilePaths.append(filePath)
            else:
                continue
    
    def returnLineWithTheBestHit(self, filePath):
        ecLineListAfterThreashold = []
        with open(filePath) as ecFile:
            line = ecFile.readlines()[15]
            if line == "\n":
                return ecLineListAfterThreashold
            splitedLine = line.split()
            if "e-" not in  splitedLine[0]:
                return ecLineListAfterThreashold
            e_valueSplited = splitedLine[0].split("-")
            if int(e_valueSplited[1]) >= self.threshold:
                ecLineListAfterThreashold.append(line.strip())
        return ecLineListAfterThreashold
    
    def thresholdingMethod(self, filePath):
        """Method which id adding files after thresholding to a list and returns it

        Args:
            filePath (str): path of a file to checked

        Returns:
            list: List of lines after thresholding
        """
        ecLineListAfterThreashold = []
        with open(filePath) as ecFile:
            for line in ecFile.readlines()[15:]:
                if line == "\n":
                    # check if 16 line is empty if yes break loop
                    break
                splitedLine = line.split()
                if "e-" not in splitedLine[0]:
                    # if e- (patern) not in line continue
                    continue
                e_valueSplited = splitedLine[0].split("-")
                if int(e_valueSplited[1]) >= self.threshold:
                    # thresholding condition
                    ecLineListAfterThreashold.append(line.strip())
                else:
                    break
        return ecLineListAfterThreashold
    
    def addIDwithECsToDictBasedOnLines(self, lines:list, ecNum):
        for searchLine in lines:
            # e_value = searchLine.split()[0]
            splitedSearchLine = searchLine.split("|")
            iD = splitedSearchLine[-2]
            if iD not in self.peptideIDecListDict:
                self.peptideIDecListDict[iD] = [ecNum]
            elif iD in self.peptideIDecListDict and ecNum not in self.peptideIDecListDict[iD]:
                self.peptideIDecListDict[iD].append(ecNum)
            

    def addToDict(self, oneHitBool):
        """Method adds KO-EC number (key) and list of peptide IDs after thresholding (value) to a dict
        """
        for filePath in self.listOfFilePaths:
            splitedPath = filePath.split(":")
            ecNum = splitedPath[-1][:splitedPath[-1].index("]")]
            if oneHitBool:
                ecLineListAfterThreashold = self.returnLineWithTheBestHit(filePath)
            else:
                ecLineListAfterThreashold = self.thresholdingMethod(filePath)
            if len(ecLineListAfterThreashold) == 0:
                continue
            else:
                self.addIDwithECsToDictBasedOnLines(ecLineListAfterThreashold, ecNum)
    
    def saveEcKoToFile(self):
        with open(self.fileName, "w") as saveFile:
            saveFile.write("Query ID\tPredicted EC number\n")
            for peptideID in self.peptideIDecListDict:
                listToStringOfECs = " ".join(self.peptideIDecListDict[peptideID])
                # listToStringOfECs = list(self.peptideIDecListDict[peptideID])
                saveFile.write(f"{peptideID}\t{listToStringOfECs}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Threshold EC number based on E-value")
    parser.add_argument("sdir", metavar="search_dir", type=str, 
                        help="Path to direcotry wih EC number hits against targeted sequence")
    parser.add_argument("-o", metavar="output", dest= "output", type=str, default="ECnumsAfterTsh.txt",
                        help="change a name of a file to save")
    parser.add_argument("-th", metavar="threshold_vale", dest="th", type=int, default=50,
                        help="""change a thresholing value, int number represents an negative index of a power of 10.
                        So -th 50 stands for E-value * 10^-50""")
    parser.add_argument("-bh", action="store_true",
                        help="Return only one the best hit from all hits. Default False -> not one hit")
    args = parser.parse_args()
    sdir = args.sdir
    output = args.output
    th = args.th
    bestHit = args.bh
    
    threshold = Threshold(sdir, output, th, bestHit)
    threshold.saveEcKoToFile()