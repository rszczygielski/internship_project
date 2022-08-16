import argparse
import numpy as np
import os
from phylodm import PhyloDM
from sklearn.metrics import silhouette_score
import pandas as pd

class SilhouetteIndex():
    def __init__(self, treeDir, clusterDir):
        self.treeDir = treeDir
        self.clusterDir = clusterDir
        self.treeFileList = []
        self.clusterFileList = []
        self.addFilesToList(".tree", self.treeDir, self.treeFileList)
        self.addFilesToList(".tsv", self.clusterDir, self.clusterFileList)

    def addFilesToList(self, ext, direcotryPath, listOfFlies):
        """Method adds all fasta format files from directory to a list
        """
        for fileInDir in os.listdir(direcotryPath):
            if fileInDir.endswith(ext):
                filePath = os.path.join(direcotryPath, fileInDir)
                listOfFlies.append(filePath)
            else:
                continue
    
    def createDictBasedOnPathList(self, pathList):
        fileNameFilePathDict = {}
        for path in pathList:
            splitedPath = path.split("/")
            fileNameWithExtension = splitedPath[-1]
            fileNameWithoutExtension = fileNameWithExtension[:fileNameWithExtension.index("]")+1]
            fileNameFilePathDict[fileNameWithoutExtension] = path
        return fileNameFilePathDict
    
    def returnFileNameBasedOnPath(self, path):
        splitedPath = path.split("/")
        fileNameWithExtension = splitedPath[-1]
        fileNameWithoutExtension = fileNameWithExtension[:fileNameWithExtension.index("]")+1]
        return fileNameWithoutExtension

    def returnDistanceMatrix(self, newickFile):
        return newickFile.dm(norm=False)
    
    def returnListOfLabels(self, newickFile):
        labels = newickFile.taxa()
        for elem in enumerate(labels):
            newElem = elem[1].strip()
            newElem = newElem[newElem.index("_")+1:]
            labels[elem[0]] = newElem
        return labels

    def returnSumOfAllSequences(self, allValuesCount):
        return allValuesCount.sum()
    
    def returnListOfPeptidesInOneCluster(self, dataSeries, peptideID):
        # return list(dataSeries.where(dataSeries.values == peptideID).dropna().index)[1:]
        return list(dataSeries.where(dataSeries.values == peptideID).dropna().index)

    def returnListOfClusterIDsToTake(self, dataSeries, allValuesCount):
        oneClusterBool = True
        listOfClusterIDsTaken = []
        biggestCluster = allValuesCount.idxmax()
        sumOfAllSequences = self.returnSumOfAllSequences(allValuesCount)
        if allValuesCount[biggestCluster] > sumOfAllSequences / 2:
            listOfClusterIDsTaken.extend(self.returnListOfPeptidesInOneCluster(dataSeries, biggestCluster))
        else:
            oneClusterBool = False
            for cluster in allValuesCount.index[:2]:
                listOfClusterIDsTaken.extend(self.returnListOfPeptidesInOneCluster(dataSeries, cluster))
        return listOfClusterIDsTaken, oneClusterBool
    
    def returnListOfClusterIDsWhichAreNotTaken(self, dataSeries, allValuesCount, oneClusterBool):
        listOfClusterIDsNotTaken = []
        if oneClusterBool:
            for cluster in allValuesCount.index[1:]:
                listOfClusterIDsNotTaken.extend(self.returnListOfPeptidesInOneCluster(dataSeries, cluster))
        else:
            for cluster in allValuesCount.index[2:]:
                listOfClusterIDsNotTaken.extend(self.returnListOfPeptidesInOneCluster(dataSeries, cluster))
        return listOfClusterIDsNotTaken
    
    def replaceDotWithUnderscore(self, peptideID):
        return peptideID.replace(".", "_")
    
    def returnDictWithPeptideIDsAndClusterLabels(self, listOfClusterIDsTaken, listOfClusterIDsNotTaken):
        peptideIDsClusterLabelsDict = {}
        for peptideID in listOfClusterIDsTaken:
            if "." in peptideID:
                peptideID = self.replaceDotWithUnderscore(peptideID)
            peptideIDsClusterLabelsDict[peptideID] = 1
        for peptideID in listOfClusterIDsNotTaken:
            if "." in peptideID:
                peptideID = self.replaceDotWithUnderscore(peptideID)
            peptideIDsClusterLabelsDict[peptideID] = 0
        return peptideIDsClusterLabelsDict

    def returnReorderdListOfLabels(self, peptideIDsClusterLabelsDict, newickFile):
        listOfReorderdClusterLabels = []
        listOfIndexesToDelete = []
        labelsInNewickFile = self.returnListOfLabels(newickFile)
        for peptideIDtuple in enumerate(labelsInNewickFile):
            if peptideIDtuple[1] not in peptideIDsClusterLabelsDict:
                listOfIndexesToDelete.append(peptideIDtuple[0])
                continue
            listOfReorderdClusterLabels.append(peptideIDsClusterLabelsDict[peptideIDtuple[1]])
        return listOfReorderdClusterLabels, listOfIndexesToDelete
     
    def returnSilhouetteIndex(self, newickFilePath, clusterFilePathDict):
        fileName = self.returnFileNameBasedOnPath(newickFilePath)
        newickFileLoaded = PhyloDM.load_from_newick_path(newickFilePath)
        dataSeries = pd.read_table(clusterFilePathDict[fileName], delim_whitespace=True, usecols=[0,1],
        index_col=[1], header=None).squeeze("columns")
        allValuesCount = dataSeries.value_counts()
        listOfClusterIDsTaken, oneClusterBool = self.returnListOfClusterIDsToTake(dataSeries, allValuesCount)
        listOfClusterIDsNotTaken = self.returnListOfClusterIDsWhichAreNotTaken(dataSeries, allValuesCount, oneClusterBool)
        peptideIDsClusterLabelsDict = self.returnDictWithPeptideIDsAndClusterLabels(listOfClusterIDsTaken, listOfClusterIDsNotTaken)
        listOfReorderdClusterLabels, listOfIndexesToDelete = self.returnReorderdListOfLabels(peptideIDsClusterLabelsDict, newickFileLoaded)
        if len(set(listOfReorderdClusterLabels)) == 1:
            return
        distanceMatrix = self.returnDistanceMatrix(newickFileLoaded)
        distanceMatrix = np.delete(distanceMatrix, listOfIndexesToDelete, axis=0)
        distanceMatrix = np.delete(distanceMatrix, listOfIndexesToDelete, axis=1)
        return round(silhouette_score(distanceMatrix, listOfReorderdClusterLabels), 2)
    
    def saveAllSilhouetteIndexes(self):
        clusterFilePathDict = self.createDictBasedOnPathList(self.clusterFileList)
        with open("silhouetteIndex.txt", "w") as silhouetteFile:
            silhouetteFile.write("KO_EC_NUMBER\tSilhouette_Index\n")
            for newickFilePath in self.treeFileList:
                silhouetteScore = self.returnSilhouetteIndex(newickFilePath, clusterFilePathDict)
                if silhouetteScore == None:
                    continue
                else:
                    fileName = self.returnFileNameBasedOnPath(newickFilePath)
                    silhouetteFile.write(f"{fileName}\t{silhouetteScore}\n")
            

if __name__ == "__main__":
    silhouetteIndex = SilhouetteIndex("/home/radeksz/Documents/AllECnumsTrees/allECnumsTrees/","/home/radeksz/bioinf/essential_files/allECnumsClustered/allMMseq2Clusters/")
    silhouetteIndex.saveAllSilhouetteIndexes()
    # print(silhouetteIndex.createDictBasedOnPathList(silhouetteIndex.clusterFileList))

    



