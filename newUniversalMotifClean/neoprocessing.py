import numpy as np
import copy   

class MatrixMakerMultiple():
    
    def __init__(self, data, mutations_to_HLA, valid_letters, pseudo_count=1, random_model=None):
        self.data = data
        self.seq_length = len(data[0][0])
        self.count_matrix = {}
        self.final_matrix = {}
        self.mutat  = mutations_to_HLA
        self.valid_letters = valid_letters
        self.letters_to_index = {x: i for i, x in enumerate(valid_letters)}
        self.pseudo = pseudo_count
        if random_model:
            self.random_model = random_model
        else:
            self.random_model = [1.0/len(self.valid_letters) for i in range(len(self.valid_letters))]
        self._build_matrix_squeleton()

    def _reset_matrix(self):
        self.count_matrix = {}
        self.final_matrix = {}
    
    def invert_string(self, string):
        return string[::-1]

    def _build_matrix_squeleton(self):
        for i in range(self.seq_length):
            self.count_matrix[i] = {}
            for mutseq in self.mutat[i]:
                self.count_matrix[i][mutseq] = {}
                for letter in self.valid_letters:
                    self.count_matrix[i][mutseq][letter] = self.pseudo
                self.count_matrix[i][mutseq]["total"] = 0
        
    def _equivals(self, position, hla):
        mutation = None
        for mutant_seq in self.mutat[position]:
            if hla in self.mutat[position][mutant_seq]:
                mutation = mutant_seq
                break
        if mutation:
            return mutation
        else:
            print("ERROR: HLA %s does not have sequence on position %s\nTALK TO THE DEVELOPER" % (hla, position))
            exit(1)
 
    def _count(self):
        for info in self.data:
            sequence_list = list(info[0])
            hla = info[1]
            for i,letter in enumerate(sequence_list):
                mut_seq = self._equivals(i,hla)
                self.count_matrix[i][mut_seq][letter] += 1
                self.count_matrix[i][mut_seq]["total"] += 1

    def build(self):
        self._count()
        self.final_matrix = copy.deepcopy(self.count_matrix)
        for num in self.final_matrix:
            for mut in self.final_matrix[num]:
                total = self.final_matrix[num][mut]["total"]  + self.pseudo * len(self.valid_letters)
                total = total * 1.0  # cross compatibility
                #debug = []
                for letter in self.valid_letters:
                    self.final_matrix[num][mut][letter] /= total
                    #debug.append(self.final_matrix[num][mut][letter])
                    self.final_matrix[num][mut][letter] /= self.random_model[self.letters_to_index[letter]]
                    self.final_matrix[num][mut][letter] = np.log2(self.final_matrix[num][mut][letter])
                    self.final_matrix[num][mut][letter] *= -1
                #print(sum(debug))

    def write_matrix(self, out, unique=True, count=False):
        if unique:
            self._parse_matrix(self.final_matrix, out)
        if count:
            self._parse_matrix(self.count_matrix, "count_%s" % out)

    def _parse_matrix(self, matrix, out):
        with open(out, "w") as inn:
            string = "Position"
            for i in range(self.seq_length):
                for mut in matrix[i]:
                    string = string + "\t%s" % str((i+1))
            string = string + "\n"
            inn.write(string)
            string = "count"
            for i in range(self.seq_length):
                for mut in matrix[i]:
                    string = string + "\t%s" % matrix[i][mut]["total"]
            string = string + "\n"
            inn.write(string)
            string = "seq"
            for i in range(self.seq_length):
                for mut in matrix[i]:
                    string = string + "\t%s" % mut
            string = string + "\n"
            inn.write(string)
            for letter in self.valid_letters:
                string = letter
                for i in range(self.seq_length):
                    for mut in matrix[i]:
                        string = string + "\t" + str(round(matrix[i][mut][letter], 3) * -1)
                string = string + "\n"
                inn.write(string)

    def puntuate_peptide(self, seq, hla):
        peptide_score = 0.0
        for position, letter in enumerate(seq):
            peptide_score += self.final_matrix[position][self._equivals(position, hla)][letter]
        return peptide_score

    def recompute_matrix(self,penalti):
        changes = 0
        for i,info in enumerate(self.data):
            hla = info[1]
            seq = info[0]
            reverse_seq = self.invert_string(seq)
            sense_score = self.puntuate_peptide(seq, hla)
            anti_sense_score = self.puntuate_peptide(reverse_seq, hla)
            if (anti_sense_score- penalti) < sense_score:
                self.data[i] = (reverse_seq, hla)
                changes += 1
        print("%s inversions of %s peptides with penalti %s" % (changes,len(self.data),penalti))
        if changes:
            self._reset_matrix()
            self._build_matrix_squeleton()
            self.build()
        return changes

    def iterative_peptide_flipping(self, n=0, penalti = 1.0):
        print("Starting iteration %s" % n)
        new_penalti = penalti - (penalti/10.0) * n
        if new_penalti < 0:
            new_penalti = 0
        changes = self.recompute_matrix(new_penalti)
        if n == 20:
            print("Maximum recursive reached, no convergence found")
            return
        if changes:
            self.iterative_peptide_flipping(n+1, penalti)
        
    def filter_outliers(self, filter_value=90, criteria="percentatge"):
        scores = []
        valid_methods = ["percentatge", "threshold"]
        if criteria == "percentatge":
            for info in self.data:
                hla = info[1]
                seq = info[0]
                score = self.puntuate_peptide(seq, hla)
                scores.append(score)
            scores = sorted(scores)
            index = int(round(len(scores)*filter_value)/100.0)
            filter_value = scores[index]
        elif criteria not in valid_methods:
            print("Error method %s not in valid_methods" % criteria)
            print("Valid filtering methods:")
            print(valid_methods)
            return
        index_to_remove = []
        for i,info in enumerate(self.data):
            hla = info[1]
            seq = info[0]
            score = self.puntuate_peptide(seq, hla)
            if score > filter_value:
                index_to_remove.append(i)
        for j in index_to_remove:
            self.data = self.data[:j] + self.data[j+1 :]
        self._reset_matrix()
        self._build_matrix_squeleton()
        self.build()
        


if __name__ == "__main__":
    import hlaizer
    import pickle
    valid = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
    d2, data = hlaizer.main()
    Motif = MatrixMakerMultiple(data, d2, valid)
    Motif.build()
    # Motif.write_matrix("DEVELOPER.txt", unique=True, count=False)
    #Motif.iterative_peptide_flipping(penalti=3)
    #Motif.write_matrix("DEVELOPER_inverted.txt", unique=True, count=False)
    #with open("inverted_motif_3.pkl","wb") as inn:
    #    pickle.dump(Motif.final_matrix, inn)
    # exemple filter outliers
    Motif.filter_outliers(filter_value=90, criteria="percentatge") # per quedaere amb el 90%
    with open("clean_motif_90_2.pkl","wb") as inn:
        pickle.dump(Motif.final_matrix, inn)
    #Motif.filter_outliers(filter_value=-5, criteria="threshold")   # per filtrar aquells que siguin mes positius que -5
                          
            
    

   








