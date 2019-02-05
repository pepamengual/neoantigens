from __future__ import print_function
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
plt.style.use("ggplot")


def readFile(name_file, transform=True):
    with open(name_file) as f:
        sequences_list = []
        dG_list = []
        for line in f:
            values = line.rstrip().split()
            seq_name = values[0]
            ic50 = values[1]
            sequences_list.append(seq_name)
            if transform:
                dG_list.append(R*T*np.log(float(ic50)/1e9))
            else:
                dG_list.append(float(ic50))
    return sequences_list, dG_list


states = ["I", "O"]
polarity = {'A':'A','R':'P','N':'P','D':'P','C':'A','E':'P','Q':'P','G':'A','H':'P','I':'A','L':'A','K':'P','M':'A','F':'A','P':'A','S':'P','T':'P','W':'A','Y':'A','V':'A'}
alphabet = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
# alphabet = ["A", "P"]
stateOrder = {state: i for i, state in enumerate(states)}
alphabetOrder = {char: i for i, char in enumerate(alphabet)}

emission = np.zeros((len(states), len(alphabet)))
transition = np.zeros((len(states), len(states)))
path = ["I", "O", "I", "O", "I", "O", "I", "O", "I"]
path = path[1:]+path[1:2]
path = ["O", "I", "I", "O", "O", "O", "I", "O", "I"]
pseudocounts = 0.0001


for i, state in enumerate(path[:-1]):
    transition[stateOrder[state], stateOrder[path[i+1]]] += 1

R = 1.9872036e-3
T = 300
sequences, dG = readFile("9-lenght-IC50-nM.txt")


for sequence in sequences:
    for state, letter in zip(path, sequence):
        # letter = polarity[letter]
        emission[stateOrder[state], alphabetOrder[letter]] += 1

for row in transition:
    tot = row.sum()
    if tot > 0:
        add = (pseudocounts*tot/(1-pseudocounts*tot))
    else:
        add = pseudocounts
    row += add

for row in emission:
    tot = row.sum()
    if tot > 0:
        add = (pseudocounts*tot/(1-pseudocounts*tot))
    else:
        add = pseudocounts
    row += add

emission /= emission.sum(axis=1).reshape((emission.shape[0], 1))
transition /= transition.sum(axis=1).reshape((transition.shape[0], 1))
emission = np.log(emission)
transition = np.log(transition)

# sequences_test, dG = readFile("sequences_file.txt", transform=False)
sequences_test = sequences

print("Peptide\tProbability\tExperimental dG(kcal/mol)")
predictions = []
for seq, dG_val in zip(sequences_test, dG):
    prob = 0.0
    for state, symbol in zip(path, seq):
        # symbol = polarity[symbol]
        prob += emission[stateOrder[state], alphabetOrder[symbol]]
    print(seq, np.exp(prob), dG_val)
    predictions.append(np.exp(prob))

slope, intercept, r_value, p_value, std_err = stats.linregress(predictions, dG)
predictions = np.array(predictions)
dG = np.array(dG)
sequences_test = np.array(sequences_test)
with open("good_probabilities.txt", "w") as fw:
    for seq in sequences_test[predictions > 1e-10]:
        fw.write("%s\n" % seq)
# plt.semilogy(dG[dG < -6], predictions[dG < -6], 'x')
# plt.semilogy(dG, predictions, 'x')
# for i, label in enumerate(sequences_test):
#     plt.text(predictions[i], dG[i], label)
# plt.gca().invert_yaxis()
# plt.ylabel("Predicted probabilities")
# plt.xlabel("Experimental (kcal/mol)")
# plt.title("R^2: %.3f" % r_value**2)
# plt.show()
