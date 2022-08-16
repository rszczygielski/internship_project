import matplotlib.pyplot as plt
import os

class Make_plot():
    def __init__(self, dirWithResults):
        self.fileList = []
        self.precisionResults = {}
        self.recallResults = {}
        self.specificityResults = {}
        self.jaccard_IndexResults = {}
        self.addFilesPathToList(".txt", dirWithResults)
        self.createDictsBasedOnFiles()
    
    def addFilesPathToList(self, ext, direcotryPath):
        """Method adds all fasta format files from directory to a list
        """
        for fileInDir in os.listdir(direcotryPath):
            if fileInDir.endswith(ext):
                filePath = os.path.join(direcotryPath, fileInDir)
                self.fileList.append(filePath)
            else:
                continue
        self.fileList.sort()

    def addValesToDicts(self, txtFile):
        splitedTxtFile = txtFile.split("^-")
        power = splitedTxtFile[1][:splitedTxtFile[1].index(".")]
        wholePower = f"10^-{power}"
        with open(txtFile) as fileToRead:
            for line in fileToRead.readlines():
                if "OVERALL PRECISION SCORE:" in line:
                    splitedLine = line.split()
                    score = float(splitedLine[-1])
                    self.precisionResults[wholePower] = score * 100
                if "OVERALL RECALL SCORE:" in line:
                    splitedLine = line.split()
                    score = float(splitedLine[-1])
                    self.recallResults[wholePower] = score * 100
                # if "OVERALL Specificity SCORE:" in line:
                #     splitedLine = line.split()
                #     score = float(splitedLine[-1])
                    self.specificityResults[wholePower] = score * 100
                if "Jaccard Index:" in line:
                    splitedLine = line.split()
                    score = float(splitedLine[-1])
                    self.jaccard_IndexResults[wholePower] = score * 100
                

    def createDictsBasedOnFiles(self):
        for txtFile in self.fileList:
            self.addValesToDicts(txtFile)
    
    def plotPrecision(self):
        x = list(self.precisionResults.keys())
        y = list(self.precisionResults.values())
        plt.plot(x,y)
        plt.xlabel("E-values scientific notation")
        plt.ylabel("Percentage (%)")
        plt.show()

    def plotRecall(self):
        x = list(self.recallResults.keys())
        y = list(self.recallResults.values())
        plt.plot(x,y)
        plt.xlabel("E-values scientific notation")
        plt.ylabel("Percentage (%)")
        plt.show()
    
    def plotSpecificity(self):
        x = list(self.specificityResults.keys())
        y = list(self.specificityResults.values())
        plt.plot(x,y)
        plt.xlabel("E-values scientific notation")
        plt.ylabel("Percentage (%)")
        plt.show()

    def plotJaccard_Index(self):
        x = list(self.jaccard_IndexResults.keys())
        y = list(self.jaccard_IndexResults.values())
        plt.plot(x,y)
        plt.xlabel("E-values scientific notation")
        plt.ylabel("Percentage (%)")
        plt.show()
    
    def plotAll(self):
        font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 30}
        # plt.rcParams.update({'font.size': 4})
        plt.style.use("fivethirtyeight")
        plt.figure(figsize=(21, 11), dpi=150)
        plt.rc('xtick', labelsize=20) 
        plt.rc('ytick', labelsize=20) 
        x = list(self.precisionResults.keys())
        plt.plot(x, list(self.precisionResults.values()),"b", label="Precision")
        plt.plot(x, list(self.recallResults.values()), "g", label="Recall")
        # plt.plot(x, list(self.specificityResults.values()), "y", label="Specificity")
        plt.plot(x, list(self.jaccard_IndexResults.values()), "r", label="Jaccard Index")
        plt.xlabel("E-values scientific notation", fontdict=font)
        plt.ylabel("Percentage (%)", fontdict=font)
        plt.ylim(0, 100)
        plt.rcParams.update({'font.size': 25})
        # plt.title("The best EC hit")
        plt.title("Prediction of a EC nums best hit after clustering", fontdict=font)
        plt.legend()
        plt.savefig("Prediction_of_a_EC_nums_best_hit_after_clustering")
        # plt.show()

if __name__ == "__main__":
    make_plot = Make_plot("/home/radeksz/bioinf/essential_files/Classification/ECnumsAfterClustering/new/OneHit/classification/")
    # make_plot.plotJaccard_Index()
    # make_plot.plotPrecision()
    # make_plot.plotRecall()
    # make_plot.plotAll()
    # print(make_plot.precisionResults)
    # print(make_plot.recallResults)
    print(make_plot.jaccard_IndexResults)


