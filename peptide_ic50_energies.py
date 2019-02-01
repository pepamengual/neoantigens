import os
import numpy as np

inputFile = "HLA-A-0201_assays.csv"
outputFile = "9-lenght-IC50-nM.txt"
foldxIE = "/home/pepamengual/neoantigenos/big_benchmark/foldxIE"

polarity = {'A':'0.83','R':'0.83','N':'0.09','D':'0.64','C':'1.48','E':'0.65','Q':'0','G':'0.1','H':'1.1','I':'3.07','L':'2.52','K':'1.6','M':'1.4','F':'2.75','P':'2.7','S':'0.14','T':'0.54','W':'0.31','Y':'2.97','V':'1.79'}
heavyatoms = {'A':'6','R':'12','N':'9','D':'9','C':'7','E':'10','Q':'10','G':'5','H':'11','I':'9','L':'9','K':'10','M':'9','F':'12','P':'8','S':'7','T':'8','W':'15','Y':'13','V':'9'}

R = 1.9872036e-3
T = 300

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


with open(outputFile, 'r') as f:
    for line in f:
        peptide = line.split(" ")[0]
        peptide = peptide.rstrip()
        polarityPeptide = 0.0
        heavyatomsPeptide = 0
        for aa in peptide:
            polarityPeptide += float(polarity[aa])
        for aa in peptide:
            heavyatomsPeptide += int(heavyatoms[aa])
        polarityPeptide = round(polarityPeptide, 3)
        ic50 = line.split(" ")[1]
        ic50 = float(ic50[0])
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
            results = "{} {} {} {} {} {} {} {} {}".format(peptide, qualitativeValue, heavyatomsPeptide, polarityPeptide, G, energy, GperHydro, GperHydroN, GperN)
            l.append(results)

with open("energies_{}".format(outputFile), "w") as f:
    f.write("Peptide qualitativeValue HeavyatomsPeptide Polarity ΔGexp FoldxBindingEnergy EnergyCorrPolarity EnergyCorrPolarityN EnergyCorrN\n")
    for k in l:
        f.write(k + '\n')


