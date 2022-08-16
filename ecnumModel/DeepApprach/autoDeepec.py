import os
import argparse

class Deepec():
    def __init__(self, directoryWithEC, newDirectoryName="allDeepec"):
        self.directoryWithEC = directoryWithEC
        self.fileList = []
        self.newDirectoryName = newDirectoryName
        self.createNewDirectory()
        self.addFilesToList()

    def addFilesToList(self):
        """Method adds all fasta format files from directory to a list
        """
        ext = ('.fasta')
        for fastaFile in os.listdir(self.directoryWithEC):
            if fastaFile.endswith(ext):
                filePath = os.path.join(self.directoryWithEC, fastaFile)
                self.fileList.append(filePath)
            else:
                continue
    
    def createNewDirectory(self):
        """Method creates a new directory to save all multiple alignments
        """
        pathOfNewDir = os.path.join(self.directoryWithEC, self.newDirectoryName)
        if not os.path.isdir(pathOfNewDir):
            os.makedirs(pathOfNewDir)
    
    def clusterFasta(self, cpu):
        """Method uses Mafft to create multuple alignments
        """
        for fastaFilePath in self.fileList:
            splitedPath = fastaFilePath.split("/")
            fastaFile = splitedPath[-1]
            bracketIndex = fastaFile.index("]")
            deepecFile = f"{fastaFile[:bracketIndex+1]}"
            newDir = os.path.join(self.newDirectoryName, deepecFile)
            pathWithoutFile = splitedPath[:-1]
            pathWithoutFile = "/".join(pathWithoutFile)
            savePath = os.path.join(pathWithoutFile, newDir)
            os.system(f"python3 /home/radeksz/bioinf/deepec/deepec.py -i {fastaFilePath} -o {savePath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Program automatically creates multiple alignments for every fasta file in directory')
    parser.add_argument("ec_dir", metavar="ec_dir", type=str,
                        help='Path to directory with fasta files to be aligned')
    parser.add_argument("-c", metavar="cpu",dest="cpu", type=int, default=8, 
                        help='set the number of parallel worker threads to <n> -> default: 8')
    args = parser.parse_args()
    ec_dir = args.ec_dir
    cpu = args.cpu

    cluster = Deepec(ec_dir)
    cluster.clusterFasta(cpu)
