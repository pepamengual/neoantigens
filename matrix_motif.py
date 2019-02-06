import numpy as np
import sys

def readMotif(d):
    with open("hla-a-0201-motif.txt", "r") as f:
        next(f)
        for line in f:
            line = line.split()
            residue = line[0]
            scores = line[1:]
            if not residue in d:
                d[residue] = scores
    return d

def readPeptides(d):
    l = []
    with open("energies_9-lenght-IC50-nM.txt", "r") as f:
        next(f)
        for line in f:
            line = line.split()
            peptide = line[0]
            qualitativeValue = line[1]
            G = float(line[2])
            Foldx = float(line[3])
            Pydock = float(line[4])
            peptideValue = 0
            for position, residue in enumerate(peptide):
                value = round(float(d[residue][position]), 2)
                peptideValue += value
            peptideValue = round(peptideValue, 2)
            if peptideValue == 0:
                continue
            FoldxValue = round(-Foldx * peptideValue, 3)
            PydockValue = round(-Pydock * peptideValue, 3)
            string = "{} {} {} {} {} {} {} {}".format(peptide, qualitativeValue, G, Foldx, Pydock, peptideValue, FoldxValue, PydockValue)
            l.append(string)

    with open("9-lenght-IC50-nM-matrix.txt", "w") as f:
        f.write("Peptide qualitativeValue ΔGexp FoldxBindingEnergy PydockBindingEnergy peptideValue FoldxValue PydockValue\n")
        for k in l:
            f.write(k + '\n')

def main():
    d = {}
    readMotif(d)
    readPeptides(d)

if __name__ == "__main__":
    main()
