import matplotlib.pyplot as plt


r'$\mathrm{10^{-20}}$'

precissionOneHitAfterClustering = {r'$\mathrm{10^{-10}}$': 49.0, r'$\mathrm{10^{-20}}$': 54.0, r'$\mathrm{10^{-30}}$': 58, r'$\mathrm{10^{-40}}$': 63.0, r'$\mathrm{10^{-50}}$': 67.0, r'$\mathrm{10^{-60}}$': 71.0, r'$\mathrm{10^{-70}}$': 74.0, r'$\mathrm{10^{-80}}$': 77.0, r'$\mathrm{10^{-90}}$': 81.0, r'$\mathrm{10^{-100}}$': 83.0, r'$\mathrm{10^{-120}}$': 85.0, r'$\mathrm{10^{-140}}$': 87.0, r'$\mathrm{10^{-150}}$': 88.0}
precissionListAfterClustering = {'10^-10': 17.0, '10^-20': 23.0, '10^-30': 28, '10^-40': 36.0, '10^-50': 43.0, '10^-60': 50.0, '10^-70': 56, '10^-80': 61.0, '10^-90': 66.0, '10^-100': 71.0, '10^-120': 78.0, '10^-140': 83.0, '10^-150': 84.0}
precissonOneHitBeforeClustering = {'10^-10': 46.0, '10^-20': 51.0, '10^-30': 55.0, '10^-40': 59.0, '10^-50': 62.0, '10^-60': 66.0, '10^-70': 70.0, '10^-80': 73.0, '10^-90': 75.0, '10^-100': 77.0, '10^-120': 81.0, '10^-140': 83.0, '10^-150': 83.0}
precissonListBeforeClustering = {'10^-10': 14.0, '10^-20': 19.0, '10^-30': 23.0, '10^-40': 28.0, '10^-50': 34.0, '10^-60': 39.0, '10^-70': 45.0, '10^-80': 51.0, '10^-90': 56.0, '10^-100': 60.0, '10^-120': 67.0, '10^-140': 73.0, '10^-150': 75.0}

recallOneHitAfterClustering = {'10^-10': 87.0, '10^-20': 88.0, '10^-30': 89.0, '10^-40': 89.0, '10^-50': 90.0, '10^-60': 90.0, '10^-70': 90.0, '10^-80': 91.0, '10^-90': 91.0, '10^-100': 92.0, '10^-120': 92.0, '10^-140': 92.0, '10^-150': 93.0}
recallListAfterClustering = {'10^-10': 86.0, '10^-20': 87.0, '10^-30': 87.0, '10^-40': 87.0, '10^-50': 88.0, '10^-60': 89.0, '10^-70': 89.0, '10^-80': 90.0, '10^-90': 90.0, '10^-100': 90.0, '10^-120': 90.0, '10^-140': 91.0, '10^-150': 91.0}
recallOneHitBeforeClustering = {'10^-10': 89.0, '10^-20': 90.0, '10^-30': 90.0, '10^-40': 90.0, '10^-50': 90.0, '10^-60': 91.0, '10^-70': 92.0, '10^-80': 92.0, '10^-90': 92.0, '10^-100': 92.0, '10^-120': 92.0, '10^-140': 92.0, '10^-150': 92.0}
recallListBeforeClustering = {'10^-10': 88.0, '10^-20': 89.0, '10^-30': 89.0, '10^-40': 89.0, '10^-50': 89.0, '10^-60': 89.0, '10^-70': 89.0, '10^-80': 90.0, '10^-90': 90.0, '10^-100': 90.0, '10^-120': 90.0, '10^-140': 90.0, '10^-150': 89.0}


jaccardIndexOneHitAfterClustering = {'10^-10': 41.2, '10^-20': 45.8, '10^-30': 49.2, '10^-40': 52.8, '10^-50': 56.5, '10^-60': 58.8, '10^-70': 60.8, '10^-80': 61.9, '10^-90': 62.3, '10^-100': 62.3, '10^-120': 57.5, '10^-140': 51.8, '10^-150': 48.3}
jaccardIndexListAfterClustering = {'10^-10': 44.3, '10^-20': 49.2, '10^-30': 52.1, '10^-40': 55.4, '10^-50': 58.8, '10^-60': 61.0, '10^-70': 62.9, '10^-80': 63.8, '10^-90': 64.1, '10^-100': 64.1, '10^-120': 59.2, '10^-140': 53.4, '10^-150': 49.9}
jaccardIndexOneHitBeforeClustering = {'10^-10': 39.0, '10^-20': 43.8, '10^-30': 47.5, '10^-40': 50.8, '10^-50': 53.6, '10^-60': 56.2, '10^-70': 58.0, '10^-80': 59.4, '10^-90': 58.7, '10^-100': 57.9, '10^-120': 52.0, '10^-140': 43.8, '10^-150': 39.3}
jaccardIndexListBeforeClustering = {'10^-10': 40.8, '10^-20': 45.7, '10^-30': 49.6, '10^-40': 52.8, '10^-50': 55.7, '10^-60': 58.0, '10^-70': 59.7, '10^-80': 60.8, '10^-90': 60.1, '10^-100': 59.2, '10^-120': 53.2, '10^-140': 44.8, '10^-150': 40.5}

deepEcPrecision = {'10^-10': 89.0, '10^-20': 89.0, '10^-30': 89.0, '10^-40': 89.0, '10^-50': 89.0, '10^-60': 89.0, '10^-70': 89.0, '10^-80': 89.0, '10^-90': 89.0, '10^-100': 89.0, '10^-120': 89.0, '10^-140': 89.0, '10^-150': 89.0}
deepEcRecall = {'10^-10': 87.0, '10^-20': 87.0, '10^-30': 87.0, '10^-40': 87.0, '10^-50': 87.0, '10^-60': 87.0, '10^-70': 87.0, '10^-80': 87.0, '10^-90': 87.0, '10^-100': 87.0, '10^-120': 87.0, '10^-140': 87.0, '10^-150': 87.0}
deepEcJaccardIndex = {'10^-10': 72.0, '10^-20': 72.0, '10^-30': 72.0, '10^-40': 72.0, '10^-50': 72.0, '10^-60': 72.0, '10^-70': 72.0, '10^-80': 72.0, '10^-90': 72.0, '10^-100': 72.0, '10^-120': 72.0, '10^-140': 72.0, '10^-150': 72.0}


combinedApprachPrecision = {'10^-10': 85.0, '10^-20': 85.0, '10^-30': 85.0, '10^-40': 85.0, '10^-50': 85.0, '10^-60': 85.0, '10^-70': 85.0, '10^-80': 85.0, '10^-90': 85.0, '10^-100': 85.0, '10^-120': 85.0, '10^-140': 85.0, '10^-150': 85.0}
combinedApprachRecall = {'10^-10': 90.0, '10^-20': 90.0, '10^-30': 90.0, '10^-40': 90.0, '10^-50': 90.0, '10^-60': 90.0, '10^-70': 90.0, '10^-80': 90.0, '10^-90': 90.0, '10^-100': 90.0, '10^-120': 90.0, '10^-140': 90.0, '10^-150': 90.0}
combinedApprachJaccardIndex = {'10^-10': 74.6, '10^-20': 74.6, '10^-30': 74.6, '10^-40': 74.6, '10^-50': 74.6, '10^-60': 74.6, '10^-70': 74.6, '10^-80': 74.6, '10^-90': 74.6, '10^-100': 74.6, '10^-120': 74.6, '10^-140': 74.6, '10^-150': 74.6}

################################################################################################
########################################PRECISION###############################################
################################################################################################
font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 30}
lineWidth = 4
# plt.rcParams.update({'font.size': 4})
print(plt.style.available)
plt.style.use("seaborn-white")
plt.figure(figsize=(21, 12), dpi=150)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)
x = precissionOneHitAfterClustering.keys()
plt.plot(x, precissionOneHitAfterClustering.values(), "b", label="One Hit After Clustering", linewidth=lineWidth)
plt.plot(x, precissionListAfterClustering.values(), "g", label="List Of Hits After Clustering", linewidth=lineWidth)
plt.plot(x, precissonOneHitBeforeClustering.values(), "r", label="One Hit Before Clustering", linewidth=lineWidth)
plt.plot(x, precissonListBeforeClustering.values(), "orange", label="List Of Hits Before Clustering", linewidth=lineWidth)

plt.plot(x, deepEcPrecision.values(), "purple", label="DeepEC", linewidth=lineWidth)
plt.plot(x, combinedApprachPrecision.values(), "turquoise", label="Combined Approach", linewidth=lineWidth)



plt.xlabel("E-value", fontdict=font)
plt.ylabel("Precision score(%)", fontdict=font)
plt.ylim(0, 100)
plt.rcParams.update({'font.size': 25})
# plt.title("The best EC hit")
# plt.title("Prediction of a EC numbers precision scores", fontdict=font)
plt.legend(loc = 4)
# plt.show()
plt.savefig("Prediction_of_a_EC_numbers_precision_scores")




################################################################################################
########################################JACCARD INDEX###########################################
################################################################################################




# plt.rcParams.update({'font.size': 4})
plt.style.use("seaborn-white")
plt.figure(figsize=(21, 12), dpi=150)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)
x = precissionOneHitAfterClustering.keys()
plt.plot(x, jaccardIndexOneHitAfterClustering.values(), "b", label="One Hit After Clustering", linewidth=lineWidth)
plt.plot(x, jaccardIndexListAfterClustering.values(), "g", label="List Of Hits After Clustering", linewidth=lineWidth)
plt.plot(x, jaccardIndexOneHitBeforeClustering.values(), "r", label="One Hit Before Clustering", linewidth=lineWidth)
plt.plot(x, jaccardIndexListBeforeClustering.values(), "orange", label="List Of Hits Before Clustering", linewidth=lineWidth)

plt.plot(x, deepEcJaccardIndex.values(), "purple", label="DeepEC", linewidth=lineWidth)
plt.plot(x, combinedApprachJaccardIndex.values(), "turquoise", label="Combined Approach", linewidth=lineWidth)

plt.xlabel("E-value", fontdict=font)
plt.ylabel("Jaccard Index score(%)", fontdict=font)
plt.ylim(0, 100)
plt.rcParams.update({'font.size': 25})
# plt.title("The best EC hit")
# plt.title("Prediction of a EC numbers Jaccard Index scores", fontdict=font)
plt.legend(loc = 4)
# plt.show()
plt.savefig("Prediction_of_a_EC_numbers_Jaccard_Index_scores")



################################################################################################
######################################RECALL#############################################
################################################################################################


# plt.rcParams.update({'font.size': 4})
plt.style.use("seaborn-white")
plt.figure(figsize=(21, 12), dpi=150)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)
x = precissionOneHitAfterClustering.keys()
plt.plot(x, recallOneHitAfterClustering.values(), "b", label="One Hit After Clustering", linewidth=lineWidth)
plt.plot(x, recallListAfterClustering.values(), "g", label="List Of Hits After Clustering", linewidth=lineWidth)
plt.plot(x, recallOneHitBeforeClustering.values(), "r", label="One Hit Before Clustering", linewidth=lineWidth)
plt.plot(x, recallListBeforeClustering.values(), "orange", label="List Of Hits Before Clustering", linewidth=lineWidth)

plt.plot(x, deepEcRecall.values(), "purple", label="DeepEC", linewidth=lineWidth)
plt.plot(x, combinedApprachRecall.values(), "turquoise", label="Combined Approach", linewidth=lineWidth)

plt.xlabel("E-value", fontdict=font)
plt.ylabel("Recall score(%)", fontdict=font)
plt.ylim(75, 100)
plt.rcParams.update({'font.size': 25})
# plt.title("The best EC hit")
# plt.title("Prediction of a EC numbers recall scores", fontdict=font)
plt.legend(loc = 4)
# plt.show()
plt.savefig("Prediction_of_a_EC_numbers_recall_scores")






# plt.plot(x, list(self.precisionResults.values()),"b", label="Precision")
# plt.plot(x, list(self.recallResults.values()), "g", label="Recall")


