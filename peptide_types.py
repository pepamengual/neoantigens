import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
import statsmodels.api as sm
from pyseqlogo.pyseqlogo import draw_logo, setup_axis
plt.style.use("ggplot")

polarity = {'A':'A','R':'P','N':'P','D':'P','C':'A','E':'P','Q':'P','G':'A','H':'P','I':'A','L':'A','K':'P','M':'A','F':'A','P':'A','S':'P','T':'P','W':'A','Y':'A','V':'A'}

inputFile = "HLA-A-0201_assays.csv"

d = {}
with open(inputFile, 'r') as f:
    for line in f:
        h_line = line.split(",")
        reference = h_line[0]
        peptide = h_line[1]
        method = h_line[2]
        assay = h_line[3]
        units = h_line[4]
        qualitativeMeasure = h_line[5]
        quantitativeValue = h_line[6]
        alleleName = h_line[7].rstrip()
        if len(peptide) == 9:
            if "X" not in peptide:
                if qualitativeMeasure:
                    if not qualitativeMeasure in d:
                        d[qualitativeMeasure] = [peptide] 
                    else:
                        d[qualitativeMeasure].append(peptide)


order = ["Negative", "Positive-Low", "Positive", "Positive-Intermediate", "Positive-High"]

index = 0
colors = ['blue','cyan','green','orange','red']

for num, category in enumerate(order):
    value = d[category]
    A = [0] * 9
    P = [0] * 9
    for peptide in value:
        for pos, residue in enumerate(peptide):
            classification = polarity[residue]
            if classification == "A":
                A[pos] += 1
            elif classification == "P":
                P[pos] += 1
    
    totalA = sum(A[0:len(A)])
    totalP = sum(P[0:len(P)])

    prop = [round((x*1.0)/y, 3) for x, y in zip(A, P)]
    print(prop, category)
    plt.plot(prop, label=category, linewidth=3, c=colors[num])

plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8, 9])
plt.ylabel('Hydrophobicity')
plt.xlabel('Peptide position')
plt.legend()
plt.savefig("image.png")
plt.show()

    

