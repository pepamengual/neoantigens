import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

inputFile = "HLA-A-0201_assays.csv"
outputFile = "9-lenght-IC50-nM.txt"
energiesFile = "energies_{}".format(outputFile)
foldxIE = "/home/pepamengual/neoantigenos/big_benchmark/foldxIE"
threshold = 15

polarity = {'A':'0.83','R':'0.83','N':'0.09','D':'0.64','C':'1.48','E':'0.65','Q':'0','G':'0.1','H':'1.1','I':'3.07','L':'2.52','K':'1.6','M':'1.4','F':'2.75','P':'2.7','S':'0.14','T':'0.54','W':'0.31','Y':'2.97','V':'1.79'}
heavyatoms = {'A':'6','R':'12','N':'9','D':'9','C':'7','E':'10','Q':'10','G':'5','H':'11','I':'9','L':'9','K':'10','M':'9','F':'12','P':'8','S':'7','T':'8','W':'15','Y':'13','V':'9'}
#order = ["Negative", "Positive-Low", "Positive", "Positive-Intermediate", "Positive-High"]
#colors = ['blue','cyan','green','orange','red']
order = ["Negative", "Positive-Low", "Positive-Intermediate", "Positive-High"]
colors = ['blue','cyan','orange','red']
#order = ["Positive"]
#colors = ["green"]

polarityClassifier = {'A':'A','R':'P','N':'P','D':'P','C':'A','E':'P','Q':'P','G':'A','H':'P','I':'A','L':'A','K':'P','M':'A','F':'A','P':'A','S':'P','T':'P','W':'A','Y':'A','V':'A'}

R = 1.9872036e-3
T = 300

def DBfilter():
    d = {}
    l = []
    with open(inputFile, 'r') as f:
        for line in f:
            h_line = line.split(",")
            referemce = h_line[0]
            peptide = h_line[1]
            method = h_line[2]
            assay = h_line[3]
            units = h_line[4]
            qualitativeValue = h_line[5]
            quantitativeValue = h_line[6]
            allele = h_line[7].rstrip()
            if len(peptide) == 9:
                if quantitativeValue:
                    if assay == "half maximal inhibitory concentration (IC50)":
                        if units == "nM":
                            if not peptide in d:
                                content = "{} {}".format(quantitativeValue, qualitativeValue)
                                d[peptide] = content
                            else:
                                del d[peptide]
    with open(outputFile, 'w') as f:
        for key, value in d.items():
            f.write("{} {}".format(key, value) + '\n')

def energyCalculator():
    with open(outputFile, 'r') as f:
        for line in f:
            peptide = line.split(" ")[0]
            peptide = peptide.rstrip()
            polarityPeptide = 0.0
            heavyatomsPeptide = 0
            polarityKeyResidues = float(polarity[peptide[1]]) + float(polarity[peptide[8]])
            for aa in peptide:
                polarityPeptide += float(polarity[aa])
            for aa in peptide:
                heavyatomsPeptide += int(heavyatoms[aa])
            polarityPeptide = round(polarityPeptide, 3)
            ic50 = line.split(" ")[1]
            ic50 = float(ic50)
            ic50 /= 1e9
            qualitativeValue = line.split(" ")[2].rstrip('\n')
            locationEnergy = "{}/Summary_{}_AC.fxout".format(foldxIE, peptide)
            with open(locationEnergy, 'r') as f2:
                last_line = f2.readlines()[-1]
                energy = last_line.split("\t")[5]
                energy = round(float(energy), 3)
                G = R * T * (np.log(ic50))
                G = round(G, 3)
                GperHydro = energy * polarityPeptide
                GperHydro = round(GperHydro, 3)
                GperHydroN = (energy * polarityPeptide) / heavyatomsPeptide
                GperHydroN = round(GperHydroN, 3)
                GperN = (energy / heavyatomsPeptide)
                GperN = round(GperN, 3)
                pKey = energy * polarityKeyResidues
                pKey = round(pKey, 3)
                pKeyN = (energy * polarityKeyResidues) / heavyatomsPeptide
                pKeyN = round(pKeyN, 3)
                results = "{} {} {} {} {} {} {} {} {} {} {}".format(peptide, qualitativeValue, heavyatomsPeptide, polarityPeptide, G, energy, GperHydro, GperHydroN, GperN, pKey, pKeyN)
                l.append(results)

    with open(energiesFile, "w") as f:
        f.write("Peptide qualitativeValue HeavyatomsPeptide Polarity ΔGexp FoldxBindingEnergy EnergyCorrPolarity EnergyCorrPolarityN EnergyCorrN\n")
        for k in l:
            f.write(k + '\n')

di = {"x": [], "y": [], "p": []} #FoldX Binding energies lower than -15Kcal/mol
di2 = {"x": [], "y": [], "p": []} #FoldX Binding energies higher than -15Kcal/mol
di3 = {"x": [], "y": [], "p": []} #FoldX Binding energies lower than -15Kcal/mol AND Apolar at position 2 or 9
di4 = {"x": [], "y": [], "p": []} #FoldX Binding energies lower than -15Kcal/mol AND any Apolar at position 2 and 9
di5 = {"x": [], "y": [], "p": []} #Double polar
di6 = {"x": [], "y": [], "p": []} #Double apolar
di7 = {"x": [], "y": [], "p": []} #Some is polar

ed = {} #energydiscarted
ea = {} #energyaccepted
dp = {} #doublepolar
da = {} #doubleapolar
sp = {} #somepolar

countdpNoHigh = 0
countdpHigh = 0
countdpNegative = 0
countdpPositive = 0

for num, name in enumerate(order):
    with open(energiesFile, "r") as f:
        next(f)
        for line in f:
            h_line = line.split(" ")
            peptide = h_line[0]
            qualitativeValue = h_line[1]
            experimentalEnergy = float(h_line[4])
            foldXBindingEnergy = float(h_line[5])
            if qualitativeValue == name:
                if foldXBindingEnergy <= -threshold:
                    if not qualitativeValue in ea:
                        ea[qualitativeValue] = 1
                    else:
                        ea[qualitativeValue] += 1
                    di["p"].append(qualitativeValue)
                    di["x"].append(experimentalEnergy)
                    di["y"].append(foldXBindingEnergy)
                    if polarityClassifier[peptide[1]] == "P":
                        di3["p"].append(qualitativeValue)
                        di3["x"].append(experimentalEnergy)
                        di3["y"].append(foldXBindingEnergy)
                        if polarityClassifier[peptide[8]] == "P": # double polar
                            di5["p"].append(qualitativeValue)
                            di5["x"].append(experimentalEnergy)
                            di5["y"].append(foldXBindingEnergy)
                            if not qualitativeValue in dp:
                                dp[qualitativeValue] = 1
                            else:
                                dp[qualitativeValue] += 1
                    if polarityClassifier[peptide[1]] == "A":
                        di4["p"].append(qualitativeValue)
                        di4["x"].append(experimentalEnergy)
                        di4["y"].append(foldXBindingEnergy)
                        if polarityClassifier[peptide[8]] == "A": #double apolar
                            di6["p"].append(qualitativeValue)
                            di6["x"].append(experimentalEnergy)
                            di6["y"].append(foldXBindingEnergy)
                            if not qualitativeValue in da:
                                da[qualitativeValue] = 1
                            else:
                                da[qualitativeValue] += 1
                            if qualitativeValue == "Positive-High":
                                countdpHigh += 1
                            else:
                                countdpNoHigh += 1
                            if qualitativeValue == "Negative":
                                countdpNegative += 1
                            if "Positive" in qualitativeValue:
                                countdpPositive += 1
                    if (polarityClassifier[peptide[1]] == "P") or (polarityClassifier[peptide[8]] == "P"): #some polar
                        di7["p"].append(qualitativeValue)
                        di7["x"].append(experimentalEnergy)
                        di7["y"].append(foldXBindingEnergy)
                        if not qualitativeValue in sp:
                            sp[qualitativeValue] = 1
                        else:
                            sp[qualitativeValue] += 1
                else:
                    di2["p"].append(qualitativeValue)
                    di2["x"].append(experimentalEnergy)
                    di2["y"].append(foldXBindingEnergy)
                    if not qualitativeValue in ed:
                        ed[qualitativeValue] = 1
                    else:
                        ed[qualitativeValue] += 1

print("Threshold: ", threshold)
print("Energy Accepted: ", ea)
print("Energy Discarted: ", ed)
print("Double Polar: ", dp)
print("Double Apolar: ", da)
print("Some Polar: ", sp)
print("Double Polar High Ratio: ", round(countdpHigh/countdpNoHigh, 3))
print("Positive / Negative Ratio: ", round(countdpPositive/countdpNegative, 3))


sns.scatterplot(x="x", y="y", hue="p", s=20, palette=colors, data=di6, alpha=1)
sns.scatterplot(x="x", y="y", hue="p", s=20, palette=colors, data=di2, alpha=0.1, legend=False)


plt.plot([-threshold, -2], [-threshold, -threshold], linewidth=2, color="black")
plt.xlabel("ΔGexp (Kcal/mol)")
plt.ylabel("FoldX Interaction Energy (Kcal/mol)")
plt.ylim(-30,10)
plt.xlim(-15,-2)
plt.savefig("scatter_9-lenght_exp_foldx_doubleapolar_{}.png".format(threshold))
plt.show()



