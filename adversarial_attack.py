import numpy as np
import random
from qiskit.quantum_info import Pauli


#returns a non-identity Pauli matrix, with optional restrictions
#potential_non_identity_paulis may be specified as a binary string, where 1s indicate a qubit position
#in the Pauli that can can be a single qubit Pauli, and 0s indicate a qubit position where the
#corresponding single-qubit Pauli is always set to the identity
#max_non_identity_paulis only has an effect if potential_non_identity_paulis is not specified
#max_non_identity_paulis indicates the maximum amount of non-identity single qubit Paulis,
#with no restriction on location of the Paulis
#min_non_identity_paulis only has an effect if potential_non_identity_paulis is not specified
#min_non_identity_paulis indicates the minimum amount of non-identity single qubit Paulis,
#with no restriction on location of the Paulis
def specific_random_paulis(size_of_pauli, min_non_identity_paulis = 0, max_non_identity_paulis = 0, potential_non_identity_paulis = [-1]):
    #return a random Pauli with at most max_non_identity_paulis number of single qubit
    #non-identity Paulis
    if max_non_identity_paulis != 0 and potential_non_identity_paulis[0] == -1:
        non_identity_paulis = random.randint(min_non_identity_paulis, max_non_identity_paulis)
        potential_non_identity_paulis = max_non_identity_paulis - non_identity_paulis
        min_identity_paulis = size_of_pauli - max_non_identity_paulis
        potential_non_identity_paulis = np.array([0] * min_identity_paulis + [1] * potential_non_identity_paulis + [2] * non_identity_paulis)
        np.random.shuffle(potential_non_identity_paulis)

    #just return a random Pauli, with no restrictions other than it can not be the identity
    if max_non_identity_paulis == 0 and potential_non_identity_paulis[0] == -1:
        potential_non_identity_paulis = [1] * size_of_pauli
    
    indexed_paulis = ["I", "X", "Y", "Z"]
    #makes sure the identity matrix is never returned
    is_identity_pauli = True
    while is_identity_pauli:
        pauli_string = ""
        if potential_non_identity_paulis[0] != -1:
            for i in range(size_of_pauli):
                if potential_non_identity_paulis[i] == 1:
                    pauli_index = random.randint(0, 3)
                    pauli_string = pauli_string + indexed_paulis[pauli_index]
                elif potential_non_identity_paulis[i] == 2:
                    pauli_index = random.randint(1, 3)
                    pauli_string = pauli_string + indexed_paulis[pauli_index]
                else:
                    pauli_string = pauli_string + "I"
                if pauli_string[-1] != "I":
                    is_identity_pauli = False
        else:
            for i in range(size_of_pauli):
                pauli_index = random.randint(0, 3)
                pauli_string = pauli_string + indexed_paulis[pauli_index]
                if pauli_string[-1] != "I":
                    is_identity_pauli = False
                
    return Pauli(pauli_string)

