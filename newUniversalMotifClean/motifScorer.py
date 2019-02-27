from sklearn import metrics
import csv
import pickle
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

dictionaryHLAFile = "dictionaryFile_HLA_A.txt"
scoreM = "clean_motif_90.pkl"
mhcDatabase = "/home/pepamengual/neoantigenos/gitneoantigens/neoantigens/database/mhc_ligand_full.csv"

colors = ['blue','cyan','orange','red']
orderViolin = ["Negative", "Positive-Low", "Positive-Intermediate", "Positive-High"]
orderViolin2Groups = ["Negative", "Positive"]
colorViolin = {"Negative": "blue", "Positive-Low": 'cyan', "Positive-Intermediate": 'orange', "Positive-High": 'red'}
colorViolin2Groups = {"Negative": "blue", "Positive": 'red'}
R = 1.9872036e-3
T = 300

def readingDictionary(dictionaryHLAFile):
    s = open(dictionaryHLAFile, 'r').read()
    saveDictionary = eval(s)
    return saveDictionary

def readingScoreMatrix(scoreM):
    with open(scoreM, 'rb') as f:
        scoreMatrix = pickle.load(f)
    return scoreMatrix

def readMHCDatabase(mhcDatabase):
    l = set()
    dictToSaveMHCInformation = {}
    with open(mhcDatabase, 'r') as f:
        next(f)
        next(f)
        for line in f:
            h_line = line.split('","')
            try:
                peptide = h_line[11]
                assay = h_line[80]
                units = h_line[81]
                qualitativeValue = h_line[83]
                ic50 = h_line[85]
                allele = h_line[95]
                ic50 = float(ic50)
                ic50 /= 1e9
                G = R * T * (np.log(ic50))
                G = round(G, 3)
            except:
                continue
            if assay == "half maximal inhibitory concentration (IC50)":
                if units == "nM":
                    if len(peptide) == 9 and not "X" in peptide:
                        if qualitativeValue in orderViolin:
                            dictToSaveMHCInformation.setdefault(allele, {})
                            if not peptide in dictToSaveMHCInformation[allele]:
                                dictToSaveMHCInformation[allele].setdefault(peptide, []).append((qualitativeValue, G))
                            else:
                                l.add((allele, peptide))
    for allele, peptide in l:
        dictToSaveMHCInformation[allele].pop(peptide)
    return dictToSaveMHCInformation

def scorer(saveDictionary, scoreMatrix, hla, dictToSaveMHCInformation):
    dictToPrintViolin = {"Prediction": [], "Qualitative": []}
    dictToPrintViolin2Groups = {"Prediction": [], "Qualitative": []}
    dictToMatrixScorerROC = {}
    dataToSave = []
    for allele, values in dictToSaveMHCInformation.items():
        if allele == hla:
            for peptide, values2 in values.items():
                for qualitativeValue, G in values2:
                    if (qualitativeValue == "Negative" and G > -7.3) or (qualitativeValue == "Positive-Low" and (G < -7.3 or G > -8.63)) or (qualitativeValue == "Positive-Intermediate" and (G < -8.63 or G > -10)) or (qualitativeValue == "Positive-High" and G < -10):
                        peptideScore = 0.0
                        for peptidePosition, peptideResidue in enumerate(peptide):
                            hlaSeqPos = saveDictionary[peptidePosition][hla]
                            scorePosition = scoreMatrix[peptidePosition][hlaSeqPos][peptideResidue]
                            peptideScore += scorePosition
                        peptideScoreRounded = round(peptideScore, 3)
                        dictToPrintViolin["Prediction"].append(peptideScoreRounded)
                        dictToPrintViolin["Qualitative"].append(qualitativeValue)
                        if qualitativeValue == "Negative" or qualitativeValue == "Positive-Low":
                            dictToPrintViolin2Groups["Prediction"].append(peptideScoreRounded)
                            dictToPrintViolin2Groups["Qualitative"].append("Negative")
                        elif qualitativeValue == "Positive-Intermediate" or qualitativeValue == "Positive-High":
                            dictToPrintViolin2Groups["Prediction"].append(peptideScoreRounded)
                            dictToPrintViolin2Groups["Qualitative"].append("Positive")
                        save = (peptide, qualitativeValue, peptideScoreRounded)
                        dataToSave.append(save)
                        dictToMatrixScorerROC.setdefault(hla, {}).setdefault(qualitativeValue, []).append(peptideScoreRounded)
    dataToSaveF = sorted(dataToSave, key=lambda tup: tup[2])
    return dictToPrintViolin, dictToPrintViolin2Groups, dataToSaveF, dictToMatrixScorerROC


def dataScorer(dictToMatrixScorerROC):
    ROCToPlot = {}
    for hla, hlaValues in dictToMatrixScorerROC.items():
        hlaScoreToROC = {"TPR": [], "FPR": []}
        for i in np.arange(20.0, -20.0, -0.1):
            TotalPositive = 0
            TotalNegative = 0
            TruePositive = 0
            FalsePositive = 0
            for qualitativeValue, mlScoreList in hlaValues.items():
                for mlScore in mlScoreList:
                    if qualitativeValue == "Positive-High" or qualitativeValue == "Positive-Intermediate":
                        TotalPositive += 1
                        if mlScore < i:
                            TruePositive += 1
                    elif qualitativeValue == "Positive-Low" or qualitativeValue == "Negative":
                        TotalNegative += 1
                        if mlScore < i:
                            FalsePositive += 1
            if TotalPositive != 0:
                TPR = round(TruePositive / TotalPositive, 3)
                hlaScoreToROC["TPR"].append(TPR)
            else:
                TPR = 0
                hlaScoreToROC["TPR"].append(TPR)
            if TotalNegative != 0:
                FPR = round(FalsePositive / TotalNegative, 3)
                hlaScoreToROC["FPR"].append(FPR)
            else:
                FPR = 0
                hlaScoreToROC["FPR"].append(FPR)
        ROCToPlot.setdefault(hla, hlaScoreToROC)
    return ROCToPlot

def ROCPlotter(ROCToPlot):
    for hla, hlaValues in ROCToPlot.items():
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(hla)
        auc = metrics.auc(hlaValues["FPR"], hlaValues["TPR"])
        aucRound = round(auc, 2)
        plt.plot(hlaValues["FPR"], hlaValues["TPR"], label="UniversalScore {}".format(aucRound))
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.legend()
        #plt.savefig("imagesROC/{}".format(hla))
        plt.show()
        plt.clf()
        plt.cla()
        plt.close()


def plotterViolin(printThis, hla):
    plt.xlabel("Peptide Score")
    plt.ylabel("Qualitative Value")
    sns.violinplot(x="Prediction", y="Qualitative", s=20, palette=colorViolin, order=orderViolin, data=printThis)
    plt.title("{}".format(hla))
    plt.savefig("images2_clean/{}_violin.png".format(hla))
    #plt.show()
    plt.clf()
    plt.cla()
    plt.close()

def plotterViolin2Groups(printThis, hla):
    plt.xlabel("Peptide Score")
    plt.ylabel("Qualitative Value")
    sns.violinplot(x="Prediction", y="Qualitative", palette=colorViolin2Groups, order=orderViolin2Groups, data=printThis)
    plt.title("{}".format(hla))
    plt.savefig("images2_clean/{}_violin2Groups".format(hla))
    #plt.show()
    plt.clf()
    plt.cla()
    plt.close()

def main():
    saveDictionary = readingDictionary(dictionaryHLAFile)
    scoreMatrix = readingScoreMatrix(scoreM)
    dictToSaveMHCInformation = readMHCDatabase(mhcDatabase)
    hla_list = []
    with open("HLA-A.pfam", 'r') as f:
        for line in f:
            h_line = line.split()
            hla = h_line[0]
            hla_list.append(hla)
    for eachHLA in hla_list:
        hla = eachHLA.split("_")[0] 
        dictToPrintViolin, dictToPrintViolin2Groups, dataToSaveF, dictToMatrixScorerROC = scorer(saveDictionary, scoreMatrix, hla, dictToSaveMHCInformation)
        with open("predictions_from_violin/{}_pred.txt".format(hla), 'w') as f:
            for line in dataToSaveF:
                line = str(line)
                f.write(line + '\n')
        ROCToPlot = dataScorer(dictToMatrixScorerROC)
        ROCPlotter(ROCToPlot)        
        plotterViolin(dictToPrintViolin, hla)
        plotterViolin2Groups(dictToPrintViolin2Groups, hla)
main()


