import argparse
import numpy as np
from sumOfPairs import ReturnIndexSequence
from logger import Logger
import os

class TotalColumnScore():
    def __init__(self, msaDirPath:str, fileName:str):
        self.msaDirPath = msaDirPath
        self.fileName = fileName
    

    def returnTotalColumnScore(self, idMsaDict):
        totalColumnScore = 0
        listOfSeq = list(idMsaDict.values())
        firstSeqMaped = list(map(str, listOfSeq[0]))
        sequenceArray = np.array([firstSeqMaped])
        for sequence in listOfSeq[1:]:
            sequenceMapped = list(map(str, sequence))
            sequenceArray = np.vstack([sequenceArray, sequenceMapped])
        Logger.settings(show_date=False, show_file_name=False)
        Logger.INFO(f"Array created, array shape = {sequenceArray.shape}")
        sequenceArray = sequenceArray.T
        for i in range (sequenceArray.shape[0]):
            if np.all(sequenceArray[i] == sequenceArray[i][0]):
                totalColumnScore += 1
        Logger.INFO(f"Total Column Score: {totalColumnScore}")
        return totalColumnScore

    
    def returnTotalColumnScoreOfAllMSAsInDir(self):
        ext = ('.aln')
        listOfStringsToSave = []
        for msaFile in os.listdir(self.msaDirPath):
            if msaFile.endswith(ext):
                filePath = os.path.join(self.msaDirPath, msaFile)
                returnID = ReturnIndexSequence(filePath)
                idMsaDict = returnID.returnPeptideDict()
                msaFile = msaFile[:msaFile.index(ext)]
                lineToSave = f"{msaFile}\t{self.returnTotalColumnScore(idMsaDict)}\t{len(idMsaDict)}"
                listOfStringsToSave.append(lineToSave)
        return listOfStringsToSave

    def saveToFile(self):
        listOfStringsToSave = self.returnTotalColumnScoreOfAllMSAsInDir()
        with open(self.fileName, "w") as saveFile:
            saveFile.write("MSAFile\tTotalColumnScore\tNumberOfSequences\n")
            for lineToSave in listOfStringsToSave:
                saveFile.write(f"{lineToSave}\n")
                Logger.settings(show_date=False, show_file_name=False)
                Logger.INFO(f"Saved {lineToSave} to file")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Program counts sum of paris score for every MSA file in chosen directory")
    parser.add_argument("msa_dir", metavar="msa_dir", type=str, 
                        help="path to the directory with all MSAs")
    parser.add_argument("-nf", metavar="name_file", dest="nf", type=str,default="TotalColumnScore.txt",
                        help="change name of a saving file -> default: TotalColumnScore.txt")
    args = parser.parse_args()
    msa_dir = args.msa_dir
    name_file = args.nf

    totalColumnScore = TotalColumnScore(msa_dir, name_file)
    totalColumnScore.saveToFile()

    # /home/radeksz/bioinf/uni-gem-rec/allECnums/allMSAs
    # TEST PATH = /home/radeksz/bioinf/uni-gem-rec/test_files