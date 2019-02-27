mhcDatabase = "/home/pepamengual/neoantigenos/gitneoantigens/neoantigens/database/mhc_ligand_full.csv"
orderViolin = ["Negative", "Positive-Low", "Positive-Intermediate", "Positive-High"]


hla_list = []
with open("HLA-A.pfam", 'r') as f:
    for line in f:
        h_line = line.split()
        hla = h_line[0]
        hla_list.append(hla)


dictToSaveMHCInformation = {}
with open(mhcDatabase, 'r') as f:
    next(f)
    next(f)
    for line in f:
        h_line = line.split('","')
        try:
            peptide = h_line[11]
            qualitativeValue = h_line[83]
            allele = h_line[95] #.rstrip()
        except:
            continue
        if len(peptide) == 9 and not "X" in peptide:
            if allele in hla_list:
                if qualitativeValue in orderViolin:
                    if not allele in dictToSaveMHCInformation:
                        dictToSaveMHCInformation[allele] = []
                        dictToSaveMHCInformation[allele].append(peptide)
                    else:
                        dictToSaveMHCInformation[allele].append(peptide)

for key, value in dictToSaveMHCInformation.items():
    with open("qualitativePeptidesFastas/{}.fasta".format(key), 'w') as f:
        for each in value:
            f.write('>' + each + '\n')
            f.write(each + '\n')


