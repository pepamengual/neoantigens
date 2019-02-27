
def main():
    alignmentFile = "HLA-A.pfam"
    dictionaryFile = "dictionaryFile_HLA_A.txt"
    mhc_file = "/home/pepamengual/neoantigenos/gitneoantigens/neoantigens/database/mhc_ligand_full.csv"
    d2 = {0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}, 7: {}, 8: {}}
    d3 = {0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}, 7: {}, 8: {}}
    hla_list = []
    
    with open(alignmentFile, 'r') as f:
        for line in f:
            h_line = line.split()
            hla_name = h_line[0]
            sequence = h_line[1]
            if not hla_name in hla_list:
                hla_list.append(hla_name)
            pos0 = "{}{}{}{}{}{}{}{}{}{}".format(sequence[4], sequence[6], sequence[32], sequence[58], sequence[62], sequence[65], sequence[158], sequence[162], sequence[166], sequence[170])
            pos1 = "{}{}{}{}{}{}{}{}{}".format(sequence[6], sequence[8], sequence[44], sequence[62], sequence[65], sequence[66], sequence[69], sequence[98], sequence[158])
            pos2 = "{}{}{}{}{}{}{}{}".format(sequence[65], sequence[69], sequence[96], sequence[98], sequence[151], sequence[154], sequence[155], sequence[158])
            pos3 = "{}{}{}{}".format(sequence[64], sequence[65], sequence[68], sequence[69])
            pos4 = "{}".format(sequence[154])
            pos5 = "{}{}{}".format(sequence[68], sequence[69], sequence[72])
            pos6 = "{}{}{}{}{}{}{}{}{}".format(sequence[72], sequence[73], sequence[76], sequence[96], sequence[113], sequence[146], sequence[151], sequence[154], sequence[155])
            pos7 = "{}{}{}{}{}{}".format(sequence[71], sequence[72], sequence[75], sequence[76], sequence[142], sequence[146])
            pos8 = "{}{}{}{}{}{}{}{}{}{}{}{}".format(sequence[72], sequence[76], sequence[79], sequence[80], sequence[83], sequence[96], sequence[115], sequence[122], sequence[123], sequence[142], sequence[145], sequence[146])
    
            l = [pos0, pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8]
            for num, val in enumerate(l):
                d2[num].setdefault(val, []).append(hla_name)
            for num, val in enumerate(l):
                d3[num].setdefault(hla_name, val)

    with open("dictionary_tmp.txt", 'w') as f:
        d2_ = str(d2)
        f.write(d2_)

    with open(dictionaryFile, 'w') as f:
        d3 = str(d3)
        f.write(d3)

    tuple_list = []
    with open(mhc_file, 'r') as f:
        next(f)
        next(f)
        for line in f:
            h_line = line.split('","')
            try:
                reference = h_line[3]
                peptide = h_line[11]
                method = h_line[79]
                assay = h_line[80]
                units = h_line[81]
                qualitativeValue = h_line[83]
                quantitativeValue = h_line[85]
                allele = h_line[95]
            except:
                continue
            if allele in hla_list:
                if len(peptide) == 9 and not "X" in peptide:
                    if qualitativeValue == "Positive-High":
                        tuple_list.append((peptide, allele))
    tuple_list = set(tuple_list)
    tuple_list = list(tuple_list)
    return d2, tuple_list

