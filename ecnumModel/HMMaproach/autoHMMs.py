import os
import argparse

class HMM():
    def __init__(self, directoryWithMSAsPath, newDirectoryName):
        self.directoryWithMSAsPath = directoryWithMSAsPath
        self.fileList = []
        self.newDirectoryName = newDirectoryName
        self.newDirPathToSaveHMMs = os.path.join(os.getcwd(), self.newDirectoryName)
        self.createNewDirectory()
        self.addFilesPathToList(".aln", self.directoryWithMSAsPath)

    def addFilesPathToList(self, ext, direcotryPath):
        """Method adds all MSA format files from directory to a list
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
        if not os.path.isdir(self.newDirPathToSaveHMMs):
            os.makedirs(self.newDirPathToSaveHMMs)
    
    def createHMMs(self, cpu):
        """Method creates HMMs on basis of multiple alignment file and saves all of hmm files into chosen 
        directory - default allHMMs dir
        """
        for msaFilePath in self.fileList:
            splitedPath = msaFilePath.split("/")
            msaFile = splitedPath[-1]
            bracketIndex = msaFile.index("]")
            hmmFileName = f"{msaFile[:bracketIndex+1]}.hmm"
            outputPath = os.path.join(self.newDirPathToSaveHMMs, hmmFileName)
            os.system(f"hmmbuild --cpu {cpu} {outputPath} {msaFilePath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Program automatically creates Hidden Markov Models for every MSA file in directory using HMMER')
    parser.add_argument("msa_dir", metavar="msa_dir", type=str,
                        help='path to directory with aln files to be aligned')
    parser.add_argument("-nd", metavar="new_dir", type=str, dest="new_dir", default="allHMMs",
                        help='path to new directory for all HMMs to be stored in - default allHMMs dir')
    parser.add_argument("-c", metavar="cpu",dest="cpu", type=str, default=os.cpu_count(),
                        help='set the number of parallel worker threads to <n> -> default: all cpus available')
    args = parser.parse_args()
    msa_dir = args.msa_dir
    new_dir = args.new_dir
    cpu = args.cpu

    hmm = HMM(msa_dir, new_dir)
    hmm.createHMMs(cpu)