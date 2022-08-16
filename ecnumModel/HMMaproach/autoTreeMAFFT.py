import os
import argparse
import shutil

class CreateTreeMAFFT():
    def __init__(self, directoryWithEC, newDirectoryName):
        self.directoryWithEC = directoryWithEC
        self.fileList = []
        self.newDirectoryName = newDirectoryName
        self.pathOfNewDir = os.path.join(self.directoryWithEC, self.newDirectoryName)
        self.createNewDirectory()
        self.addFilesToList(".fasta", self.directoryWithEC)
        self.copyFilesToNewDir()
        self.fileList = []
        self.addFilesToList(".fasta", self.pathOfNewDir)
    
    def addFilesToList(self, ext, direcotryPath):
        """Method adds all fasta format files from directory to a list
        """
        for fileInDir in os.listdir(direcotryPath):
            if fileInDir.endswith(ext):
                filePath = os.path.join(direcotryPath, fileInDir)
                self.fileList.append(filePath)
            else:
                continue

    def copyFilesToNewDir(self):
        for filePath in self.fileList:
            shutil.copy(filePath, self.pathOfNewDir)

    def createNewDirectory(self):
        """Method creates a new directory to save all multiple alignments
        """
        if not os.path.isdir(self.pathOfNewDir):
            os.makedirs(self.pathOfNewDir)
    
    def createMSAs(self, cpu):
        """Method uses Mafft to create multuple alignments
        """
        for fastaFilePath in self.fileList:
            splitedPath = fastaFilePath.split("/")
            fastaFile = splitedPath[-1]
            bracketIndex = fastaFile.index("]")
            treeFile = f"{fastaFile[:bracketIndex+1]}"
            savePath = os.path.join(self.pathOfNewDir, treeFile)
            os.system(f"mafft  --thread {cpu} --retree 0 --treeout --globalpair --reorder {fastaFilePath} > {savePath}")
    
    def removeFilesFromDir(self, ext, direcotryPath):
        """Method removes all files with specified extension from directory
        """
        for delFile in os.listdir(direcotryPath):
            if delFile.endswith(f".{ext}"):
                filePath = os.path.join(direcotryPath, delFile)
                os.remove(filePath)
            else:
                continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Program automatically creates distance tree for every fasta file in directory using MAFFT')
    parser.add_argument("ec_dir", metavar="ec_dir", type=str,
                        help='Path to directory with fasta files to be aligned')
    parser.add_argument("-nd", metavar="new_dir", dest="new_dir", type=str, default="allECnumsTrees",
                        help="change name of a saving direcotry, default: allECnumsTrees")
    parser.add_argument("-c", metavar="cpu",dest="cpu", type=str, default=os.cpu_count(),
                        help='set the number of parallel worker threads to <n>, default: all cpus available')
    args = parser.parse_args()
    ec_dir = args.ec_dir
    new_dir = args.new_dir
    cpu = args.cpu

    msa = CreateTreeMAFFT(ec_dir, new_dir)
    msa.createMSAs(cpu)
    # msa.removeFilesFromDir("tree", msa.pathOfNewDir)
