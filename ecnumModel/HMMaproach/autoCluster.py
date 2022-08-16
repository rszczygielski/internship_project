import os
import argparse
import shutil

class Cluster():
    def __init__(self, directoryWithEC, newDirectoryName):
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
        """Method creates a new directory to save all cluster files
        """
        pathOfNewDir = os.path.join(self.directoryWithEC, self.newDirectoryName)
        if not os.path.isdir(pathOfNewDir):
            os.makedirs(pathOfNewDir)
    
    def returnSavePathBasedOnFastaFilePath(self, fastaFilePath):
        splitedPath = fastaFilePath.split("/")
        fastaFile = splitedPath[-1]
        bracketIndex = fastaFile.index("]")
        clstrFile = f"{fastaFile[:bracketIndex+1]}"
        newDir = os.path.join(self.newDirectoryName, clstrFile)
        pathWithoutFile = splitedPath[:-1]
        pathWithoutFile = "/".join(pathWithoutFile)
        return os.path.join(pathWithoutFile, newDir)
    
    def clusterFasta(self, cdHitBool, mmseq2Bool, seqIdentity, cpu):
        """Method uses CD-HIT to cluster previously created MSAs
        """
        for fastaFilePath in self.fileList:
            savePath = self.returnSavePathBasedOnFastaFilePath(fastaFilePath)
            if cdHitBool:
                os.system(f"cd-hit -i {fastaFilePath} -o {savePath} -c {seqIdentity} -n 4 -T {cpu}")
            if mmseq2Bool:
                os.system(f"mmseqs easy-cluster {fastaFilePath} {savePath} --min-seq-id {seqIdentity} tmp")
        if mmseq2Bool:
            cwd = os.getcwd()
            path = os.path.join(cwd, "tmp")
            shutil.rmtree(path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Program automatically cluster every multiple alignments file in directory (conda required)')
    parser.add_argument("ec_dir", metavar="ec_dir", type=str,
                        help='Path to directory with fasta files to be aligned')
    parser.add_argument("-cd", action="store_true",
                        help="Use CD-HIT clustering software")
    parser.add_argument("-msq", action="store_true", 
                        help="Use MMseqs2 clustering software")
    parser.add_argument("-id", metavar="seq_identity", dest="seq_identity", type=float, default=0.62,
                        help="Sequence identity threshold, type = float, default -> 0.62")
    parser.add_argument("-nd", metavar="new_dir", dest="new_dir", type=str, default="allClusters",
                        help="change name of a saving direcotry, default: allClusters")
    parser.add_argument("-c", metavar="cpu",dest="cpu", type=int, default=0, 
                        help='set the number of parallel worker threads to <n> -> default: 0 - all cpus')
    args = parser.parse_args()
    ec_dir = args.ec_dir
    cdHit = args.cd
    mmseq2 = args.msq
    seq_idt = args.seq_identity
    n_dir = args.new_dir
    cpu = args.cpu

    cluster = Cluster(ec_dir, n_dir)
    cluster.clusterFasta(cdHit, mmseq2,seq_idt, cpu)