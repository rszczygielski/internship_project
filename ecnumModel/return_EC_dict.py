from logger import Logger
import time
import os

class LineTypes():
    DEFINITION = "DEFINITION"
    GENES = "GENES"
    END = "///"
    REFERENCE = "REFERENCE"
    BRACKET = "("
    EC = "EC:"
    ENTRY = "ENTRY"

class ReturnEcIDs():
    def __init__(self, KOfile):
        self.koFile = open(KOfile)
        self.ecNumIDsDict = {}
        self.start_time = time.time()
        self.path = os.getcwd() + f"/allECnums.txt"
        self.ecNum = ""
        self.koNum = ""
        self.peptideIDslist = []
        self.idBool = False
    
    def deleteBrackets(self, idList:list):
        """Method deleates brackets from line

        Args:
            idList (list): List of peptides IDs in a line with or without brackets

        Returns:
            (list) : List of peptides IDs without brackets
        """
        listWithoutBrackets = []
        for pepID in idList:
            if LineTypes.BRACKET in pepID:
                bracketIndex = pepID.index("(")
                pepID = pepID[:bracketIndex]
                listWithoutBrackets.append(pepID)
            else:
                listWithoutBrackets.append(pepID)
        return listWithoutBrackets
    
    def returnListOfECnumsIfmanyECnumsInLine(self, ecNum:str):
        """Method returns list of EC numbers if many EC numbers in line

        Args:
            ecNum (str): many EC numbers in line

        Returns:
            (list): List of EC numbers
        """
        ecNumList = []
        ecNum = ecNum[1:-1]
        ecNums = ecNum.split()
        for elem in ecNums:
            if "-" in elem:
                continue
            if "EC:" in elem:
                ecNumList.append(f"[{elem}]")
            else:
                ecNumList.append(f"[EC:{elem}]")
        return ecNumList

    def returnECnumIfECnumLine(self, line:str):
        """Method returns EC number if EC number in line or list of EC numbers if many EC numbers in line

        Args:
            line (str): DEFINITION line in a ko file

        Returns:
            (str, list or None): EC numebr (str), List of EC numbers (list) or None if EC number not in line
        """
        if LineTypes.EC in line:
            line = line.strip()
            ecIndex = line.index("EC:")
            ecNum = line[ecIndex-1:]
            if ecNum.count(".") > 3:
                ecNumList = self.returnListOfECnumsIfmanyECnumsInLine(ecNum)
                return ecNumList
            if "-" in ecNum:
                return
            return ecNum

    def saveECnumToFile(self, withIDs=False):
        """Saves all EC numbers in a file

        Args:
            withIDs (bool, optional): Saves EC numbers with their peptides IDs in file. Defaults to False.
        """
        with open(self.path, "a") as allECnumFile:
            if withIDs:
                for key, value in self.ecNumIDsDict.items():
                    allECnumFile.write(f"{key}:\n")
                    for peptideID in value:
                        idToSave = f"\t{peptideID}"
                        allECnumFile.write(f"{idToSave} \n")
                Logger.INFO("Saved EC numbers with IDs into file")
            else:
                for key in self.ecNumIDsDict:
                    allECnumFile.write(key)
                    allECnumFile.write("\n")
                Logger.INFO("Saved EC numbers without IDs into file")
        allECnumFile.close()
    
    def addKeyValueToDictAndResetVariables(self):
        """Adds EC number (key) and list of all related peptide IDs (value) to a dictionary, and resets values:
        self.ecNum = ""
        self.peptideIDslist = []
        self.idBool = False
        if "///" or "REFERENCE" in line
        """
        if type(self.ecNum) == list:
            for elem in self.ecNum:
                koNumEcNum = f"{self.koNum}-{elem}"
                self.ecNumIDsDict[koNumEcNum] = self.peptideIDslist       
        else:
            koNumEcNum = f"{self.koNum}-{self.ecNum}"
            self.ecNumIDsDict[koNumEcNum] = self.peptideIDslist
        self.ecNum = ""
        self.koNum = ""
        self.peptideIDslist = []
        self.idBool = False
    
    def changeIdBoolIntoTrue(self, splitedLine:str):
        """Method adds IDs from first line of headline GENES and changes idBool into True - tells the program to 
        start adding protein IDs into list

        Args:
            splitedLine (str): Splited line form a file
        """
        idInFirstLine = splitedLine[2:]
        idWithoutBrackets = self.deleteBrackets(idInFirstLine)
        self.peptideIDslist.extend(idWithoutBrackets)
        self.idBool = True
    
    def returnEcIDsDict(self):
        """Returns dictionary with EC number (key) and list of all related peptide IDs (value)

        Returns:
            dictionary: Dictionary with EC number (key) and list of all related peptide IDs (value)
        """
        for line in self.koFile.readlines():
            splitedLine = line.split()
            if LineTypes.ENTRY == splitedLine[0]:
                self.koNum = splitedLine[1]
            if LineTypes.DEFINITION == splitedLine[0]:
                self.ecNum = self.returnECnumIfECnumLine(line)
                if self.ecNum == None:
                    continue
            if LineTypes.GENES == splitedLine[0] and self.ecNum != None:
                self.changeIdBoolIntoTrue(splitedLine)
                continue
            if self.idBool:
                idInLine = splitedLine[1:]
                idWithoutBrackets = self.deleteBrackets(idInLine)
                self.peptideIDslist.extend(idWithoutBrackets)
            if LineTypes.REFERENCE == splitedLine[0]:
                self.idBool = False
                continue
            if LineTypes.END == splitedLine[0] and self.ecNum != None:
                self.addKeyValueToDictAndResetVariables()
        self.koFile.close()
        Logger.settings(show_date=False, show_file_name=False)
        Logger.INFO("EC number: ID list dict was created", "This action took",
        str(round(time.time() - self.start_time, 2)), "seconds to run")        
        return self.ecNumIDsDict

if __name__ == "__main__":
    ecNum = ReturnEcIDs("/home/radeksz/bioinf/essential_files/ko")
    ecNum.returnEcIDsDict()
    ecNum.saveECnumToFile(True)