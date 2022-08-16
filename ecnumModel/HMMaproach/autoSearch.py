import os
import argparse

class HMMsearch():
    def __init__(self, directoryWithHMMsPath, targetSequenceDatabase):
        self.directoryWithHMMsPath = directoryWithHMMsPath
        self.fileList = []
        self.targetSequenceDatabase = targetSequenceDatabase
        self.pathOfNewDir = os.path.join(self.directoryWithHMMsPath,"searchDir")
        self.createNewDirectory()
        self.addFilesPathToList(".hmm", self.directoryWithHMMsPath)

    def addFilesPathToList(self, ext, direcotryPath):
        """Method adds all HMM files from directory to a list
        """
        for fileInDir in os.listdir(direcotryPath):
            if fileInDir.endswith(ext):
                filePath = os.path.join(direcotryPath, fileInDir)
                self.fileList.append(filePath)
            else:
                continue
    
    def createNewDirectory(self):
        """Method creates a new directory to save all query searches against target 
        """
        if not os.path.isdir(self.pathOfNewDir):
            os.makedirs(self.pathOfNewDir)
    
    def searchTargetGenome(self, cpu):
        """Method creates output search files on basis of hmm files and saves those files into searchDir
        """
        for hmmFilePath in self.fileList:
            splitedPath = hmmFilePath.split("/")
            hmmFile = splitedPath[-1]
            bracketIndex = hmmFile.index("]")
            searchFileName = f"{hmmFile[:bracketIndex+1]}.out"
            outputPath = os.path.join(self.pathOfNewDir, searchFileName)
            os.system(f"hmmsearch --cpu {cpu} {hmmFilePath} {self.targetSequenceDatabase} > {outputPath}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Program automatically searches profiles HMMs in direcotry against target database/sequence')
    parser.add_argument("hmm_dir", metavar="hmm_dir", type=str,
                        help='path to directory with hmm files to quary')
    parser.add_argument("t_db", metavar="t_db", type=str,
                        help='path to target data base/sequence file')
    parser.add_argument("-c", metavar="cpu",dest="cpu", type=str, default=os.cpu_count(),
                        help='set the number of parallel worker threads to <n>, default: all cpus available')
    args = parser.parse_args()
    hmm_dir =args.hmm_dir
    t_db = args.t_db
    cpu = args.cpu

    hmmsearch = HMMsearch(hmm_dir, t_db)
    hmmsearch.searchTargetGenome(cpu)

# /home/radeksz/bioinf/fasta_files/E.coli-genome.fasta
# /home/radeksz/bioinf/fasta_files/S.cerevisiae-genome.fasta