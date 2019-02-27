import hlaizer
from neoprocessing import *
valid = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
d2, data = hlaizer.main()
Motif = MatrixMakerMultiple(data, d2, valid)
Motif.build()
scoreMatrix = Motif.final_matrix
import pickle
with open('scoreMatrix.pkl', 'wb') as f:
    pickle.dump(scoreMatrix, f)
