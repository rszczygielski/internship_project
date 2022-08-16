import pandas as pd
from logger import Logger


class MetaboliteStruct():
    def __init__(self, metaboliteID, name, reference, formula, charge, mass, inChI, inChIKey, smiles):
        self.metaboliteID = metaboliteID
        self.name = name
        self.reference = reference
        self.formula = formula
        self.charge = charge
        self.mass = mass
        self.inChI = inChI
        self.inChIKey = inChIKey
        self.smiles = smiles
    
    def __str__(self):
        return f"{self.metaboliteID}\t{self.name}\t{self.reference}\t{self.formula}\t{self.charge}\t{self.mass}\t{self.inChI}\t{self.inChIKey}\t{self.smiles}"


class MapMetabolites():
    def __init__(self, mnxIDsCombines, mnxMetaboliteFile):
        self.mnxIDsCombines = mnxIDsCombines
        self.mnxMetaboliteFile = mnxMetaboliteFile
        self.listOfMetaboliteStructs = []
        Logger.settings(show_date=False, show_file_name=False)
    

    def returnMetaboliteIDsSet(self):
        dataFrame = pd.read_table(self.mnxIDsCombines, usecols=["mnx_equation"])
        listOfListEquationSplited = dataFrame["mnx_equation"].str.split().to_list()
        listOfMetabolitesExtended = []
        for row in listOfListEquationSplited:
            listOfMetabolitesExtended.extend(row)
        listOfMetabolites = [metabolite for metabolite in listOfMetabolitesExtended if "MNXM" in metabolite]
        return set(listOfMetabolites)
    
    def addOneMetaboliteStructToList(self, metaboliteID, mnxDataFrame):
        metaboliteIDWithoutBrucket = metaboliteID[:metaboliteID.index("[")]
        rowsWithEC = mnxDataFrame.loc[mnxDataFrame["ID"] == metaboliteIDWithoutBrucket]
        numpyArrayRowsListed = rowsWithEC.values.astype(str)
        for row in numpyArrayRowsListed: 
            name = row[1]
            reference = row[2] 
            formula = row[3]
            charge = row[4]
            mass = row[5]
            inChI = row[6]
            inChIKey = row[7]
            smiles = row[8]
        self.listOfMetaboliteStructs.append(MetaboliteStruct(metaboliteID, name, reference, formula, charge, mass, inChI, inChIKey, smiles))
        Logger.INFO(f"{metaboliteID} metabolite struct added to list")    
    
    def addMetabolitesToList(self):
        mnxDataFrame = pd.read_table(self.mnxMetaboliteFile, header=351).drop_duplicates()
        setOfMetaboliteIDs = self.returnMetaboliteIDsSet()
        for metaboliteID in setOfMetaboliteIDs:
            self.addOneMetaboliteStructToList(metaboliteID, mnxDataFrame)

    def saveResultsToFile(self, output):
        with open(output, "w") as resultFile:
            resultFile.write("metaboliteID\tname\treference\tformula\tcharge\tmass\tinChI\tinChIKey\tsmiles\n")
            for metaboliteStruct in self.listOfMetaboliteStructs:
                resultFile.write(f"{metaboliteStruct}\n")

if __name__ == "__main__":
    mapMetabolites = MapMetabolites("/home/radeksz/bioinf/uni-gem-rec/All_MNX_IDs_Combined.txt", "/home/radeksz/Downloads/chem_prop.tsv")
    # mapMetabolites.returnMetaboliteIDsSet()
    mapMetabolites.addMetabolitesToList()
    mapMetabolites.saveResultsToFile("Metabolites_INFO.txt")
    