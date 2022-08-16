import pandas as pd
import sys
import os
sys.path.insert(1, "/home/radeksz/bioinf/uni-gem-rec/MSAquality")
from sumOfPairs import ReturnIndexSequence

class ClusterSorter():
    def __init__(self, clusterDir, ecNumDir):
        self.clusterDir = clusterDir
        self.ecNumDir = ecNumDir
        self.numberOfBiggerClusters = 0
        self.allNumberOfClusters = 0
        self.fileList = []
        self.addFilesPathToList(ext=".tsv", direcotryPath=self.clusterDir)

    def addFilesPathToList(self, ext, direcotryPath):
        """Method adds all fasta format files from directory to a list
        """
        for fileInDir in os.listdir(direcotryPath):
            if fileInDir.endswith(ext):
                filePath = os.path.join(direcotryPath, fileInDir)
                self.fileList.append(filePath)
            else:
                continue

    def returnListOfPeptidesInOneCluster(self, dataSeries, peptideID):
        return list(dataSeries.where(dataSeries.values == peptideID).dropna().index)
    
    def returnSumOfAllSequences(self, allValuesCount):
        return allValuesCount.sum()
    
    def returnListWithBiggestClusterOrTwoTheBiggestClusters(self, allValuesCount, halfBool):
        listWithTheBiggestClusters = [allValuesCount.idxmax()]
        if not halfBool:
            listWithTheBiggestClusters.append(allValuesCount[1:].idxmax())
        return listWithTheBiggestClusters
    
    def checkIfTheBiggestClusterContainsMoreThanHalfSeq(self, allValuesCount):
        biggestCluster = allValuesCount.idxmax()
        sumOfAllSequences = self.returnSumOfAllSequences(allValuesCount)
        self.allNumberOfClusters += 1
        if allValuesCount[biggestCluster] > sumOfAllSequences / 2:
            self.numberOfBiggerClusters += 1
            return True
        else:
            return False
    
    def returnFastaECnumPathBasedOnMMMseq2ECnumCluster(self, mmseq2ClusterFile, dirWithECums):
        splitedClusterPath = mmseq2ClusterFile.split("/")
        extIndex = splitedClusterPath[-1].index("_cluster.tsv")
        clusterFile = splitedClusterPath[-1][:extIndex]
        fastaFile = f"{clusterFile}.fasta"
        return os.path.join(dirWithECums, fastaFile)
    
    def replaceEcNumsAfterClusters(self, listWithTheBiggestClusters, dataSeries, fastaFilePath, ecNumDict):
        with open(fastaFilePath, "w") as fastaFile:
            for cluster in listWithTheBiggestClusters:
                listOfPeptideIDs = self.returnListOfPeptidesInOneCluster(dataSeries, cluster)
                for peptideID in listOfPeptideIDs:
                    seq = ""
                    if peptideID in ecNumDict:
                        seq += f">{peptideID}\n"
                        seq += ecNumDict[peptideID]
                        fastaFile.write(seq)
                        fastaFile.write("\n")
    
    def returnEcNumDict(self, fastaFilePath):
        returnID = ReturnIndexSequence(fastaFilePath)
        return returnID.returnPeptideDict()

    def sortClusters(self):
        for clusterFile in self.fileList:
            dataSeries = pd.read_table(clusterFile, delim_whitespace=True, usecols=[0,1], index_col=[1], header=None).squeeze("columns")
            if len(dataSeries) == 0:
                continue
            fastaFilePath = self.returnFastaECnumPathBasedOnMMMseq2ECnumCluster(clusterFile, self.ecNumDir)
            ecNumDict = self.returnEcNumDict(fastaFilePath)
            allValuesCount = dataSeries.value_counts()
            halfBool = self.checkIfTheBiggestClusterContainsMoreThanHalfSeq(allValuesCount)
            listWithTheBiggestClusters = self.returnListWithBiggestClusterOrTwoTheBiggestClusters(allValuesCount, halfBool)
            # self.replaceEcNumsAfterClusters(listWithTheBiggestClusters, dataSeries, fastaFilePath, ecNumDict)
     

if __name__ == "__main__":
    clusterSorter = ClusterSorter("/home/radeksz/bioinf/essential_files/allECnumsCopy/allMMseq2Clusters/", "/home/radeksz/bioinf/essential_files/allECnumsCopy/")
    clusterSorter.sortClusters()
    print(f"Number of all clusters {clusterSorter.allNumberOfClusters}")
    print(f"Number of cluster with sequences more than 50% {clusterSorter.numberOfBiggerClusters}")
    print(f"Procet of cluster with sequences more than 50% {round((clusterSorter.numberOfBiggerClusters/clusterSorter.allNumberOfClusters)*100, 2)}%")
    