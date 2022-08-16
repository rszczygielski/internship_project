import blosum
import sys
import time
import os
import argparse
sys.path.insert(1, "/home/radeksz/bioinf/uni-gem-rec/ecnumModel")
from logger import Logger


class ReturnIndexSequence():
    def __init__(self, seqFile):
        self.seqFile = open(seqFile)
        self.idSequenceDict = {}

    def returnPeptideDict(self):
        """Return dictionary with peptide IDs (key) with their sequence (value)

        Returns:
            (dict) : key: peptide ID, value: sequence
        """
        peptideSequence = ""
        peptideID = ""
        for line in self.seqFile.readlines():
            if ">" in line:
                self.idSequenceDict[peptideID] = peptideSequence
                peptideSequence = ""
                peptideID = line.split()[0][1:]
            else:
                peptideSequence += line.strip()
        if "" in self.idSequenceDict:
            del self.idSequenceDict[""]
        self.idSequenceDict[peptideID] = peptideSequence
        self.seqFile.close()
        return self.idSequenceDict

class ReturnBlosumScore():
    def __init__(self, msaDirPath:str, nameOfAFile, blosumType, gapPenalty):
        self.msaDirPath = msaDirPath
        self.nameOfAFile = nameOfAFile
        self.blosumType = int(blosumType)
        self.gapPenalty = float(gapPenalty)
        self.start_time = time.time()
    
    def returnScoreOfOneMSA(self, idMsaDict):
        overalScore = 0
        matrix = blosum.BLOSUM(self.blosumType, self.gapPenalty)
        listOfValues = list(idMsaDict.values())
        for i in range(len(listOfValues)):
            sequence1 = listOfValues[i]
            for j in range(i+1, len(listOfValues)):
                sequence2 = listOfValues[j]
                for k in range(len(sequence1)):
                    overalScore += matrix[f"{sequence1[k]}{sequence2[k]}"]
        Logger.settings(show_date=False, show_file_name=False)
        Logger.INFO("Sum of pairs score created", "This action took", 
        str(round(time.time() - self.start_time, 2)), "seconds to run")
        return round(overalScore, 2)
    
    def returnScoreOfAllMSAsInDir(self):
        ext = ('.aln')
        listOfStringsToSave = []
        for msaFile in os.listdir(self.msaDirPath):
            if msaFile.endswith(ext):
                filePath = os.path.join(self.msaDirPath, msaFile)
                returnID = ReturnIndexSequence(filePath)
                idMsaDict = returnID.returnPeptideDict()
                msaFile = msaFile[:msaFile.index(ext)]
                lineToSave = f"{msaFile}\t{self.returnScoreOfOneMSA(idMsaDict)}\t{len(idMsaDict)}"
                listOfStringsToSave.append(lineToSave)
        return listOfStringsToSave

    def saveToFile(self):
        listOfStringsToSave = self.returnScoreOfAllMSAsInDir()
        with open(self.nameOfAFile, "a") as saveFile:
            saveFile.write("MSAFile\tSumOfPairs\tNumberOfSequences\n")
            for lineToSave in listOfStringsToSave:
                saveFile.write(f"{lineToSave}\n")
                Logger.settings(show_date=False, show_file_name=False)
                Logger.INFO(f"Saved {lineToSave} to file")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Program counts sum of paris score for every MSA file in chosen directory")
    parser.add_argument("msa_dir", metavar="msa_dir", type=str, 
                        help="path to the directory with all MSAs")
    parser.add_argument("-nf", metavar="name_file", dest="nf", type=str,default="SumOfPairScore.txt",
                        help="change name of a saving file -> default: SumOfPairScore.txt")
    parser.add_argument("-bt", metavar="blosum_type", dest="bt", default=62,
                        help="change blosum type -> default: BLOSUM 62")
    parser.add_argument("-gp", metavar="gap_panalty", dest="gp", default=-1.53,
                        help="change gap penalty -> default: -1.53")
    args = parser.parse_args()
    msa_dir = args.msa_dir
    name_file = args.nf
    blosum_type = args.bt
    gap_panalty = args.gp

    returnBolsumScore = ReturnBlosumScore(msa_dir, name_file, blosum_type, gap_panalty)
    returnBolsumScore.saveToFile()