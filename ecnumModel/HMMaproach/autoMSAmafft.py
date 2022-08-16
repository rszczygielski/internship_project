import os
import argparse

class MSA():
    def __init__(self, directoryWithEC, newDirectoryName):
        self.directoryWithEC = directoryWithEC
        self.fileList = []
        self.newDirectoryName = newDirectoryName
        self.createNewDirectory()
        self.addFilesPathToList(".fasta", self.directoryWithEC)

    def addFilesPathToList(self, ext, direcotryPath):
        """Method adds all fasta format files from directory to a list
        """
        for fileInDir in os.listdir(direcotryPath):
            if fileInDir.endswith(ext):
                filePath = os.path.join(direcotryPath, fileInDir)
                self.fileList.append(filePath)
            else:
                continue
    
    def createNewDirectory(self):
        """Method creates a new directory to save all multiple alignments
        """
        pathOfNewDir = os.path.join(self.directoryWithEC, self.newDirectoryName)
        if not os.path.isdir(pathOfNewDir):
            os.makedirs(pathOfNewDir)

    def createMSAs(self, cpu, bl):
        """Method uses Mafft to create multuple alignments
        """
        for fastaFilePath in self.fileList:
            splitedPath = fastaFilePath.split("/")
            fastaFile = splitedPath[-1]
            bracketIndex = fastaFile.index("]")
            msaFile = f"{fastaFile[:bracketIndex+1]}.aln"
            newDir = os.path.join(self.newDirectoryName, msaFile)
            pathWithoutFile = splitedPath[:-1]
            pathWithoutFile = "/".join(pathWithoutFile)
            savePath = os.path.join(pathWithoutFile, newDir)
            os.system(f"mafft --auto --thread {cpu} --bl {bl} --maxiterate 1000 {fastaFilePath} > {savePath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Program automatically creates multiple alignments for every fasta file in directory using MAFFT')
    parser.add_argument("ec_dir", metavar="ec_dir", type=str,
                        help='Path to directory with fasta files to be aligned')
    parser.add_argument("-nd", metavar="new_dir", dest="new_dir", type=str, default="allMSAs",
                        help="change name of a saving direcotry, default: allMSAs")
    parser.add_argument("-c", metavar="cpu",dest="cpu", type=str, default=os.cpu_count(),
                        help='set the number of parallel worker threads to <n>, default: all cpus available')
    parser.add_argument("-bl", metavar="bl", dest="bl", type=int, default=62,
                        help="BLOSUM number matrix is used. number=30, 45, 62 or 80, default: 62")
    args = parser.parse_args()
    ec_dir = args.ec_dir
    new_dir = args.new_dir
    cpu = args.cpu
    bl = args.bl

    msa = MSA(ec_dir, new_dir)
    msa.createMSAs(cpu, bl)
