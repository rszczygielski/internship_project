import os
from phylodm import PhyloDM
from sklearn.cluster import KMeans
import numpy as np 
import matplotlib.pyplot as plt
import statistics

class AutoPlot():
    def __init__(self, treeDir):
        self.treeDir = treeDir
        self.fileList = []
        self.addFilesPathToList(".tree", self.treeDir)
    
    def returnMedian(self, row):
        return statistics.median(row)

    def addFilesPathToList(self, ext, direcotryPath):
        """Method adds all files with defined exntension from directory to a list
        """
        for fileInDir in os.listdir(direcotryPath):
            if fileInDir.endswith(ext):
                filePath = os.path.join(direcotryPath, fileInDir)
                self.fileList.append(filePath)
            else:
                continue
    
    def ReturnCluserPeptideIDs(self, clusterIndeciesMatrix, labelsPeptideMatrix):
        return np.take(labelsPeptideMatrix, clusterIndeciesMatrix)
    
    def returnClusterIndicesNumpyMatrix(self, clustNum, labels_array):
        return np.where(labels_array == clustNum)[0]
    
    def returnListOfmeansForEveryRow(self, distanceMatrix):
        listOfmeansForEveryRow = []
        meanMarix = np.apply_along_axis(self.returnMedian, axis=1, arr=distanceMatrix)
        for menan in meanMarix:
            listOfmeansForEveryRow.append(menan)
        return listOfmeansForEveryRow
    
    def returnListOfClusterIDs(self, identified_clusters):
        listOfClusterIDs = []
        for iD in identified_clusters:
            listOfClusterIDs.append(iD)
        return listOfClusterIDs

    def showPlot(self, listOfClusterIDs, listOfmeansForEveryRow, kmeans, ecNumber):
        plt.style.use("fivethirtyeight")
        plt.figure(figsize=(12, 6), dpi=150)
        plt.scatter(listOfClusterIDs, listOfmeansForEveryRow)
        plt.xlabel("Cluster IDs")
        plt.ylabel("Medain for every row")
        plt.tight_layout()
        # plt.ylim(0, 1.5)
        plt.figtext(.4, .8, ecNumber)
        plt.figtext(.4, .7, f"Number clusters equal to 1: {len(self.returnClusterIndicesNumpyMatrix(1, kmeans.labels_))}")
        plt.figtext(.4, .6, f"Number clusters equal to 0: {len(self.returnClusterIndicesNumpyMatrix(0, kmeans.labels_))}")
        plt.show()
        
    def returnEcNumBasedOnPath(self, treeFile):
        splitedTreeFile = treeFile.split("-")
        closingBracketIndex = splitedTreeFile[-1].index("]")
        return splitedTreeFile[-1][:closingBracketIndex+1]
        
    def autoPlot(self):
        for treeFile in self.fileList:
            pdm = PhyloDM.load_from_newick_path(treeFile)
            distanceMatrix = pdm.dm(norm=False)
            kmeans = KMeans(n_clusters=2, n_init=20)
            identified_clusters = kmeans.fit_predict(distanceMatrix)
            listOfmeansForEveryRow = self.returnListOfmeansForEveryRow(distanceMatrix)
            listOfClusterIDs = self.returnListOfClusterIDs(identified_clusters)
            ecNumber = self.returnEcNumBasedOnPath(treeFile)
            print(f"Ploting {ecNumber}")
            self.showPlot(listOfClusterIDs, listOfmeansForEveryRow, kmeans, ecNumber)


if __name__ == "__main__":
    autoPlot = AutoPlot("/home/radeksz/bioinf/essential_files/allECnumsClustered/allClusterdECnumsTrees/")
    autoPlot.autoPlot()