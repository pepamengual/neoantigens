import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
import statsmodels.api as sm
from pyseqlogo.pyseqlogo import draw_logo, setup_axis


def peptideClassifier(inputFile, d8, d9, d10):
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
            if len(peptide) == 8:
                if "X" not in peptide:
                    if qualitativeMeasure:
                        if not qualitativeMeasure in d8:
                            d8[qualitativeMeasure] = [peptide] 
                        else:
                            d8[qualitativeMeasure].append(peptide)
            if len(peptide) == 9:
                if "X" not in peptide:
                    if qualitativeMeasure:
                        if not qualitativeMeasure in d9:
                            d9[qualitativeMeasure] = [peptide]
                        else:
                            d9[qualitativeMeasure].append(peptide)
            if len(peptide) == 10:
                if "X" not in peptide:
                    if qualitativeMeasure:
                        if not qualitativeMeasure in d10:
                            d10[qualitativeMeasure] = [peptide]
                        else:
                            d10[qualitativeMeasure].append(peptide)
    return d8, d9, d10


polarity = {'A':'A','R':'P','N':'P','D':'P','C':'A','E':'P','Q':'P','G':'A','H':'P','I':'A','L':'A','K':'P','M':'A','F':'A','P':'A','S':'P','T':'P','W':'A','Y':'A','V':'A'}
order = ["Negative", "Positive-Low", "Positive", "Positive-Intermediate", "Positive-High"]
polarity = {'A':'A','R':'P','N':'P','D':'P','C':'A','E':'P','Q':'P','G':'A','H':'P','I':'A','L':'A','K':'P','M':'A','F':'A','P':'A','S':'P','T':'P','W':'A','Y':'A','V':'A'}
colors = ['blue','cyan','green','orange','red']

def plotter(d):
    index = 0
    for num, category in enumerate(order):
        value = d[category]
        peptideLenght = len(value[0])
        A = [0] * peptideLenght
        P = [0] * peptideLenght
        for peptide in value:
            for pos, residue in enumerate(peptide):
                classification = polarity[residue]
                if classification == "A":
                    A[pos] += 1
                elif classification == "P":
                    P[pos] += 1

        for n, i in enumerate(A):
            if i == 0:
                A[n] = 1
        for n, i in enumerate(P):
            if i == 0:
                P[n] = 1
        prop = [round((x*1.0)/y, 3) for x, y in zip(A, P)]
        print("\n")
        print("Count of peptide lenght {}:".format(peptideLenght))
        print(category, A, P, prop)
        plt.plot(prop, label=category, linewidth=3, c=colors[num])

    oldTicks = list(range(peptideLenght))
    newTicks = list(range(1, peptideLenght+1))
    plt.xticks(oldTicks, newTicks)
    plt.ylabel('Hydrophobicity')
    plt.xlabel('Peptide position')
    plt.legend()
    plt.savefig("peptide_lenght_{}.png".format(peptideLenght))
    plt.show()


def main():
    inputFile = "HLA-A-0201_assays.csv"
    d8 = {}
    d9 = {}
    d10 = {}
    peptideClassifier(inputFile, d8, d9, d10)    
    plotter(d8)
    plotter(d9)
    plotter(d10)


if __name__ == "__main__":
    main()
