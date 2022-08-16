import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

class PlotMsaQA():
    def __init__(self):
        self.font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}

    def plotTotalColumnScore(self):
        
        totalColumnScoreDataFrame = pd.read_table("/home/radeksz/bioinf/essential_files/TotalColumnScore.txt", index_col="MSAFile")
        zeroFilter = totalColumnScoreDataFrame["TotalColumnScore"] != 0
        totalColumnScoreDataFrame.where(zeroFilter, inplace=True)
        sortedScore = totalColumnScoreDataFrame.sort_values("NumberOfSequences")
        listOfTotalCoulumnScore = list(sortedScore["TotalColumnScore"])
        listOfNumberOfSequences = list(sortedScore["NumberOfSequences"])
        plt.figure(figsize=(10, 7), dpi=150)
        plt.rc('xtick', labelsize=20) 
        plt.rc('ytick', labelsize=20)
        plt.plot(listOfNumberOfSequences, listOfTotalCoulumnScore, ".")
        plt.xlabel("Number Of Sequences", fontdict=self.font)
        plt.ylabel("Total Column Score", fontdict=self.font)
        plt.ylim(0, 1000)
        plt.xlim(0, 1000)
        plt.savefig("Total_Column_Score")
        plt.show()
    
    def plotTotalEntropy(self):
        entropyDataFrame = pd.read_table("/home/radeksz/bioinf/essential_files/TotalEntropyScore.txt", index_col="MSAFile")
        sortedScore = entropyDataFrame.sort_values("NumberOfSequences")
        listOfTotalEntropy = list(sortedScore["TotalEntropyScore"])
        listOfNumberOfSequences = list(sortedScore["NumberOfSequences"])
        plt.figure(figsize=(10, 7), dpi=150)
        plt.rc('xtick', labelsize=20) 
        plt.rc('ytick', labelsize=20)
        plt.plot(listOfNumberOfSequences, listOfTotalEntropy, ".")
        plt.xlabel("Number Of Sequences", fontdict=self.font)
        plt.ylabel("Total Entropy Score", fontdict=self.font)
        plt.ylim(0, 6000)
        plt.xlim(0, 10000)
        plt.savefig("Total_Entropy")
        plt.show()

if __name__ == "__main__":
    plotMSA = PlotMsaQA()
    plotMSA.plotTotalColumnScore()
    plotMSA.plotTotalEntropy()