import sys
import os
from sklearn.cluster import KMeans
sys.path.insert(1, "/home/radeksz/bioinf/uni-gem-rec/MSAquality")
from sumOfPairs import ReturnIndexSequence
from phylodm import PhyloDM
import numpy as np

class TreeCluster():
    def __init__(self, directoryWithECpath):
        self.directoryWithECpath = directoryWithECpath
        self.fileList = []
        self.addFilesPathToList(".tree", self.directoryWithECpath)
    
    def addFilesPathToList(self, ext, direcotryPath):
        """Method adds all fasta format files from directory to a list
        """
        for fileInDir in os.listdir(direcotryPath):
            if fileInDir.endswith(ext):
                filePath = os.path.join(direcotryPath, fileInDir)
                self.fileList.append(filePath)
            else:
                continue

    def returnLabelsMatrix(self, newickFile):
        # newickFile = PhyloDM.load_from_newick_path(newickPath)
        # distanceMatrix = newickFile.dm(norm=False)
        labels = newickFile.taxa()
        for elem in enumerate(labels):
            newElem = elem[1].strip()
            labels[elem[0]] = newElem
        return np.array(labels)
    
    def returnClusterIndexesNumpyMatrix(self, clustNum, labels_array): #numpy 
        return np.where(labels_array == clustNum)[0]

    def returnCluserPeptideIDs(self, clusterIndeciesMatrix, labelsPeptideMatrix):
        return np.take(labelsPeptideMatrix, clusterIndeciesMatrix)

    def returnFastaECnumPathBasedOnNewickECnum(self, newickFile):
        extIndex = newickFile.index(".tree")
        fastaFile = newickFile[:extIndex]
        return os.path.join(self.directoryWithECpath, fastaFile)
    
    def countClusterAndReturnBiggerOne(self, identified_clusters):
        numberOfClusters = set(identified_clusters)
        dictOfelementsInClusters = {}
        # counting clusters number loop
        for cluster in numberOfClusters:
            numberOfEmenInCluster = np.count_nonzero(identified_clusters == cluster)
            dictOfelementsInClusters[cluster] = numberOfEmenInCluster
        return max(dictOfelementsInClusters, key=dictOfelementsInClusters.get)

    def returnMatrixOfPeptideIDsAfterClustering(self, newickFilePath):
        newickFile = PhyloDM.load_from_newick_path(newickFilePath)
        distanceMatrix = newickFile.dm(norm=False)
        labelsMatrix = self.returnLabelsMatrix(newickFile)
        kmeans = KMeans(n_clusters=2, n_init=50)
        identified_clusters = kmeans.fit_predict(distanceMatrix)
        biggerCluster = self.countClusterAndReturnBiggerOne(identified_clusters)
        clusterIndexesNumpyMatrix = self.returnClusterIndexesNumpyMatrix(biggerCluster, kmeans.labels_)
        return self.returnCluserPeptideIDs(clusterIndexesNumpyMatrix, labelsMatrix)
    
    def returnDictWithIDseqBasedOnNewickFile(self, newickFile):
        fastaPath = self.returnFastaECnumPathBasedOnNewickECnum(newickFile)
        returnIndexSequence = ReturnIndexSequence(fastaPath)
        return returnIndexSequence.returnPeptideDict()

    # loop through peptide IDs add peptide id to the seq and n/ and seq and save it in already existing fasta file


if __name__ == "__main__":
    treeCluster = TreeCluster("/home/radeksz/bioinf/essential_files/TenECnums/allECnumsTrees")
    np.set_printoptions(threshold=sys.maxsize)
    # print(treeCluster.returnMatrixOfPeptideIDsAfterClustering("/home/radeksz/bioinf/essential_files/TenECnums/allECnumsTrees/K00003-[EC:1.1.1.3].fasta.tree"))
    print(treeCluster.returnDictWithIDseqBasedOnNewickFile("K00003-[EC:1.1.1.3].fasta.tree"))



    
    

