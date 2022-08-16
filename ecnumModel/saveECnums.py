from return_EC_dict import ReturnEcIDs
from return_IDseq_dict import ReturnIndexSequence
from logger import Logger
import os
import time
import argparse

class ECnumSaver():
    def __init__(self):
        self.start_time = time.time()
        self.path = os.getcwd()     

    def saveToFile(self, seq:str, koEcNum:str):
        """Method saves all sequenase for specific EC number

        Args:
            seq (str): Fasta formated sequence > index: and a sequence in new line
        """
        path = self.path + f"/{koEcNum}.fasta"
        with open(path, "a") as fastaFile:
            fastaFile.write(seq)
            fastaFile.write("\n")
            fastaFile.close()
        
    def changePathToSave(self, newSavePath):
        """Method enable to chanege a saving  file path

        Args:
            path (str): Path of saving sequences for specific EC numebr in a file
        """
        self.path = newSavePath
        Logger.INFO(f"Changed saving path to {self.path}")

    def createFastaSeqAndSaveIt(self, koEcNum, peptideDict, ecNumIDsDict):
        """Method opens a file and creates a fasta formated sequence like >index: {amino sequence} and saves it into file
        """
        idList = ecNumIDsDict[koEcNum]
        for peptideID in idList:
            seq = ""
            if peptideID in peptideDict:
                seq += f">{peptideID}\n"
                seq += peptideDict[peptideID]
                self.saveToFile(seq, koEcNum)
            else:
                continue
        Logger.INFO(f"Created {koEcNum}.fasta in a current directory", "This action took", 
        str(round(time.time() - self.start_time, 2)), 
        "seconds to run")
    
    def saveAllEc(self, ecNumIDsDict):
        """Method saves all EC numbers with KEGG Orthology

        Args:
            ecNumIDsDict (dict): Dictionary with EC numbers as a key and peptide sequences as a value
        """
        for koEcNum in ecNumIDsDict:
            self.createFastaSeqAndSaveIt(koEcNum, peptideDict, ecNumIDsDict)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Program saves EC number with a KEGG Orthology and all sequences related to those EC numbers')
    parser.add_argument("ko_path", metavar="ko_path", type=str,
                        help='Path to KEGG Orthology file')
    parser.add_argument("sq_path", metavar="sq_path", type=str,
                        help="Path to file with IDs and sequences")
    parser.add_argument("-np", metavar="new_path", dest="new_path",
                        help="Change path to save from deafalt - current direcory")
    parser.add_argument("-ec", metavar="ec_num", dest="ecNum", 
                        help="Enter a single EC number and save all sequences related to this EC number in fasta file - defalt save all EC numbers")
    args = parser.parse_args()
    ko_path = args.ko_path
    sq_path = args.sq_path
    new_path = args.new_path
    ecNum = args.ecNum

    ecNumIDs = ReturnEcIDs(ko_path)
    ecNumSequence = ReturnIndexSequence(sq_path)
    peptideDict = ecNumSequence.returnPeptideDict()
    ecNumIDsDict = ecNumIDs.returnEcIDsDict()

    fastaFile = ECnumSaver()
    if new_path:
        fastaFile.changePathToSave(new_path)
    if ecNum:
        fastaFile.createFastaSeqAndSaveIt(ecNum, peptideDict, ecNumIDsDict)
    else:
        fastaFile.saveAllEc(ecNumIDsDict)
    