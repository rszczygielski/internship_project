from logger import Logger
import time

class ReturnIndexSequence():
    def __init__(self, seqFile):
        self.seqFile = seqFile
        self.idSequenceDict = {}
        self.start_time = time.time()

    def returnPeptideDict(self):
        """Return dictionary with peptide IDs (key) with their sequence (value)

        Returns:
            (dict) : key: peptide ID, value: sequence
        """
        peptideSequence = ""
        peptideID = ""
        with open(self.seqFile) as seqFile:
            for line in seqFile.readlines():
                if ">" in line:
                    self.idSequenceDict[peptideID] = peptideSequence
                    peptideSequence = ""
                    indexName = line.split()[0]
                    peptideID = indexName[indexName.index(":")+1:]
                else:
                    peptideSequence += line.strip()
        if "" in self.idSequenceDict:
            del self.idSequenceDict[""]
        self.idSequenceDict[peptideID] = peptideSequence
        Logger.settings(show_date=False, show_file_name=False)
        Logger.INFO("ID: sequence dict was created", "This action took", 
        str(round(time.time() - self.start_time, 2)), "seconds to run")
        return self.idSequenceDict

if __name__ == "__main__":
    peptideDict = ReturnIndexSequence("/home/radeksz/bioinf/essential_files/prokaryotes.pep")