import argparse
import numpy as np
from sumOfPairs import ReturnIndexSequence
from logger import Logger
import os
import math

class Entropy():
    def __init__(self, msaDirPath:str, fileName:str, threshold:float):
        self.msaDirPath = msaDirPath
        self.fileName = fileName
        self.threshold = threshold

    def returnTransposedSequenceArray(self, idMsaDict):
        listOfSeq = list(idMsaDict.values())
        firstSeqMaped = list(map(str, listOfSeq[0]))
        sequenceArray = np.array([firstSeqMaped])
        for sequence in listOfSeq[1:]:
            sequenceMapped = list(map(str, sequence))
            sequenceArray = np.vstack([sequenceArray, sequenceMapped])
        Logger.settings(show_date=False, show_file_name=False)
        Logger.INFO(f"Array created, array shape = {sequenceArray.shape}")
        # now sequenceArray is created
        return sequenceArray.T
    
    def returnEntropy(self, idMsaDict):
        """Method returns total Entropy Score of MSA, based on previously created dict of this exact MSA

        Args:
            idMsaDict (dict): Peptide ID - MSA dictionary

        Returns:
            float: Total Entorpy Score
        """
        totalEntropyScore = 0
        sequenceArrayTransposed = self.returnTransposedSequenceArray(idMsaDict)
        # chnages the orientaion of rows and columns for easier logic
        numOfColumns = sequenceArrayTransposed.shape[0]
        for i in range(numOfColumns):
            columnEntropy = 0
            # set of caracters in a column
            setOfColumn = set(sequenceArrayTransposed[i])
            # gaps frequency in a column
            gapsFreq = np.count_nonzero(sequenceArrayTransposed[i] == "-") / numOfColumns
            if gapsFreq > self.threshold:
                continue
            if "-" in setOfColumn:
                setOfColumn.remove("-")
                # deletes gaps from set
            for amino in setOfColumn:
                m = np.count_nonzero(sequenceArrayTransposed[i] == amino)
                # counts how many amino present in column sequenceArrey transverted so sequenceArrey.shape[1] == number of rows
                aminoFreq = m/sequenceArrayTransposed.shape[1]
                numerator = aminoFreq * math.log2(aminoFreq)
                entoropy =  -(numerator / (1 - gapsFreq))
                # counts an entropy
                columnEntropy += entoropy
            totalEntropyScore += columnEntropy
            # after looping through amino in column adding columnEntropy to totalEntropyScore
        totalEntropyScore = round(totalEntropyScore, 3)
        Logger.INFO(f"Total Entropy: {totalEntropyScore}")
        return totalEntropyScore
    
    def returnEntropyOfAllMSAsInDir(self):
        """Method returns a list of strings to save from every MSA file ext(aln) in chosen directory

        Returns:
            list: List of strings to save
        """
        ext = ('.aln')
        listOfStringsToSave = []
        for msaFile in os.listdir(self.msaDirPath):
            if msaFile.endswith(ext):
                filePath = os.path.join(self.msaDirPath, msaFile)
                returnID = ReturnIndexSequence(filePath)
                idMsaDict = returnID.returnPeptideDict()
                msaFile = msaFile[:msaFile.index(ext)]
                lineToSave = f"{msaFile}\t{self.returnEntropy(idMsaDict)}\t{len(idMsaDict)}"
                listOfStringsToSave.append(lineToSave)
        return listOfStringsToSave

    def saveToFile(self):
        """Method saves all previously created strings to file
        """
        listOfStringsToSave = self.returnEntropyOfAllMSAsInDir()
        with open(self.fileName, "w") as saveFile:
            saveFile.write("MSAFile\tTotalEntropyScore\tNumberOfSequences\n")
            for lineToSave in listOfStringsToSave:
                saveFile.write(f"{lineToSave}\n")
                Logger.settings(show_date=False, show_file_name=False)
                Logger.INFO(f"Saved {lineToSave} to file")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Program counts entropy score for every MSA file in chosen directory")
    parser.add_argument("msa_dir", metavar="msa_dir", type=str, 
                        help="path to the directory with all MSAs")
    parser.add_argument("-nf", metavar="name_file", dest="nf", type=str,default="TotalEntropyScore.txt",
                        help="change name of a saving file -> default: TotalEntropyScore.txt")
    parser.add_argument("-th", metavar="threshold", dest="th", type=float, default=0.8,
                        help="set a threshold for the gaps frequency -> default: 0.8")
    args = parser.parse_args()
    msa_dir = args.msa_dir
    name_file = args.nf
    threshold = args.th

    returnEntropy = Entropy(msa_dir, name_file, threshold)
    returnEntropy.saveToFile()