from qiskit import QuantumRegister
from qiskit.quantum_info import random_statevector, Statevector
from qiskit import QuantumCircuit
import numpy as np
import random
from qiskit import quantum_info
from qiskit.quantum_info import Clifford
import math
from sortedcontainers import SortedList

def my_sdg(qc, qubit):
    qc.s(qubit)
    qc.s(qubit)
    qc.s(qubit)

#swap qubits at position qubit_1 and qubit_1 + 1
def my_swap(qc, qubit_1):
    qc.cx(qubit_1 + 1,qubit_1)
    qc.h(qubit_1)
    qc.h(qubit_1 + 1)
    qc.cx(qubit_1 + 1,qubit_1)
    qc.h(qubit_1)
    qc.h(qubit_1 + 1)
    qc.cx(qubit_1 + 1,qubit_1)
    
def my_circ_simple(qc, starting_position, ending_position):
    if starting_position < ending_position:
        for i in range(starting_position,ending_position):
            qc.swap(i,i+1)
    elif starting_position > ending_position:
        for i in range(starting_position, ending_position, -1):
            qc.swap(i-1,i)
    
def my_circ(qc, starting_position, ending_position):
    if starting_position < ending_position:
        for i in range(starting_position,ending_position):
            my_swap(qc, i)
    elif starting_position > ending_position:
        for i in range(starting_position, ending_position, -1):
            my_swap(qc, i - 1)

def my_cx(qc, control, target):
    if control > target:
        my_circ(qc, target, control - 1)
        qc.cx(control, control - 1)
        my_circ(qc, control - 1, target)
    elif control < target:
        my_circ(qc, target, control + 1)
        qc.h(control)
        qc.h(control + 1)
        qc.cx(control + 1,control)
        qc.h(control)
        qc.h(control + 1)
        my_circ(qc, control + 1, target)   
        
def my_cy(qc, control, target):
    if control > target:
        my_circ(qc, target, control - 1)
        my_sdg(qc, control - 1)
        qc.cx(control, control - 1)
        qc.s(control - 1)
        my_circ(qc, control - 1, target)
    elif control < target:
        my_circ(qc, target, control + 1)
        my_sdg(qc, control + 1)
        qc.h(control)
        qc.h(control + 1)
        qc.cx(control + 1, control)
        qc.h(control)
        qc.h(control + 1)
        qc.s(control + 1)
        my_circ(qc, control + 1, target)

def my_cz(qc, control, target):
    if control > target:
        my_circ(qc, target, control - 1)
        qc.h(control - 1)
        qc.cx(control, control - 1)
        qc.h(control - 1)
        my_circ(qc, control - 1, target)
    elif control < target:
        my_circ(qc, target, control + 1)
        qc.h(control)
        qc.cx(control + 1,control)
        qc.h(control)
        my_circ(qc, control + 1, target)


def quick_n_choose_k(new_n, new_k, prev_n, prev_k, prev_n_choose_k):
    new_n_choose_k = 0
    if prev_n_choose_k == 0:
        if new_n == new_k:
            return 1, new_n, new_k
        else:
            return 0, new_n, new_k
    elif new_n - prev_n == 1 and prev_k == new_k:
        new_n_choose_k = prev_n_choose_k / (prev_n - (prev_k-1))
        new_n_choose_k = new_n_choose_k * new_n
    elif prev_n - new_n == 1 and prev_k == new_k:
        new_n_choose_k = prev_n_choose_k / prev_n
        new_n_choose_k = new_n_choose_k * (prev_n - prev_k)
    elif prev_n - new_n == 1 and prev_k - new_k == 1:
        new_n_choose_k = prev_n_choose_k / (prev_n)
        new_n_choose_k = new_n_choose_k * prev_k
    return new_n_choose_k, new_n, new_k


def combinatorial_number_system(n,k,N):
    elements_to_pick = []

    #for finding the first c_i, must start at c_i = k - 1, then increase c_i
    #continually until c_i choose i exceeds N. Then, use c_i - 1 as the position
    #of the ith element in the list
    c_i = k - 1
    i = k
    c_i_choose_i = 0
    prev_c_i_choose_i = c_i_choose_i
    while c_i_choose_i <= N:
        prev_c_i_choose_i = c_i_choose_i
        if c_i_choose_i == N:
            c_i = c_i + 1
            break
        c_i = c_i + 1
        c_i_choose_i, _, _ = quick_n_choose_k(c_i, i, c_i - 1, i, c_i_choose_i)
    c_i_choose_i = prev_c_i_choose_i
    N = N - c_i_choose_i
    c_i = c_i - 1
    elements_to_pick.append(c_i)        
    for i in range(k - 1, 0, -1):
        if N == 0:
            elements_to_pick.append(i-1)
        else:
            c_i = c_i - 1
            c_i_choose_i, _, _ = quick_n_choose_k(c_i, i, c_i + 1, i + 1, c_i_choose_i)
            while c_i_choose_i > N:
                c_i = c_i - 1
                c_i_choose_i, _, _ = quick_n_choose_k(c_i, i, c_i + 1, i, c_i_choose_i)
            N = N - c_i_choose_i
            elements_to_pick.append(c_i)

    return elements_to_pick

def determine_max_subkey_b(num_total_qubits, num_signature_qubits):
    #calculating key space in base 10, subtract 1
    #subtracting 1 because the first key is 00000.... 0
    #this is in effect, the largest 'value' the key can take on
    
    num_unique_outputs = 1
    m_over_d = int((num_total_qubits - num_signature_qubits)/num_signature_qubits)
    for i in range(num_signature_qubits):
        curr_n_choose_k_base = math.comb(num_total_qubits - i * m_over_d, m_over_d)
        num_unique_outputs = num_unique_outputs * curr_n_choose_k_base
    
    max_key_base_2_length = math.ceil(math.log2(num_unique_outputs))
    max_key_base_10 = (2 ** max_key_base_2_length) - 1
    
    return max_key_base_2_length, max_key_base_10

def validate_subkey_b_size(key, num_total_qubits, num_signature_qubits, verbose = False):
    max_key, _ = determine_max_subkey_b(num_total_qubits, num_signature_qubits)
    
    if verbose:
        if max_key != len(key):
            print("Invalid key with " + str(len(key)) + " bits. Key must have " + str(max_key) + " bits.")
        else:
            print("Valid key with " + str(len(key)) + " bits for a packet with " + str(num_total_qubits) + " total qubits and " + str(num_signature_qubits) + " signature qubits.")
    
    return max_key == len(key)

def subkey_b_mapping(num_total_qubits, num_signature_qubits, subkey_b, key_is_binary = True):
    output = [num_signature_qubits] * num_total_qubits
    unpicked_elements = SortedList(range(num_total_qubits))
    num_data_qubits = num_total_qubits - num_signature_qubits
    m_over_d = int(num_data_qubits / num_signature_qubits)
    integer = 0
    if key_is_binary == True:
        #this is done in O(n(O(nlogn) + O(n)) = O(n^2logn)
        for bit in subkey_b:
            integer = integer * 2 #done in O(nlogn)
            integer = integer + bit #done in O(n)
    else:
        integer = subkey_b

    for i in range(num_signature_qubits):
        curr_n_choose_k_base = math.comb(num_total_qubits - i * m_over_d, m_over_d)
        integer, curr_N = divmod(integer, curr_n_choose_k_base)
        curr_combination = combinatorial_number_system(num_total_qubits - i * m_over_d,m_over_d,curr_N)
        for element in curr_combination:
            output[unpicked_elements.pop(element)] = i
    
    return output

#some code taken from here https://stackoverflow.com/questions/3973685/python-homework-converting-any-base-to-any-base
def subkey_a_mapping(key,gates_to_use,num_qubits,key_is_binary = True):
    original_base = 2
    new_base = 3
    
    #the key is in base 2, with n digits.

    integer = 0
    if key_is_binary == True:
        #this is done in O(n(O(nlogn) + O(n)) = O(n^2logn)
        for bit in key:
            integer = integer * original_base #done in O(nlogn)
            integer = integer + bit #done in O(n)
    else:
        integer = key

    #integer 
    #this is done O(n^2logn)
    while integer:
        integer, gate = divmod(integer, new_base) #done in O(nlogn)
        gates_to_use.append(gate)
    
    for i in range(num_qubits - len(gates_to_use)):
        gates_to_use.append(0)
    gates_to_use.reverse()
        
    
    #gates_to_use = [0] * (num_qubits - len(gates_to_use)) + gates_to_use

def determine_max_subkey_a(num_qubits):
    #calculating key space in base 10, subtract 1
    #subtracting 1 because the first key is 00000.... 0
    #this is in effect, the largest 'value' the key can take on
    num_unique_outputs = 3**num_qubits
    
    max_key_base_2_length = math.ceil(math.log2(num_unique_outputs))
    max_key_base_10 = 2**(max_key_base_2_length) - 1
    
    return max_key_base_2_length, max_key_base_10

def validate_subkey_a_size(key, num_qubits, verbose = False):
    max_key, _ = determine_max_subkey_a(num_qubits)
    
    if verbose:
        if max_key != len(key):
            print("Invalid key with " + str(len(key)) + " bits. Key must have " + str(max_key) + " bits.")
        else:
            print("Valid key with " + str(len(key)) + " bits for a packet with " + str(num_qubits) + " qubits.")
    
    return max_key == len(key)

#Assume num_data_qubits is divisible by num_signature_qubits
def block_clifford_encryption(num_total_qubits, num_signature_qubits,subkey_a=None,subkey_b=None,subkey_a_base=2):
    subkey_a_gates = []
    if subkey_a == None:
        for i in range(num_total_qubits):
            subkey_a_gates.append(random.randint(0, 2))
    else:
        if subkey_a_base == 3:
            subkey_a_gates = subkey_a
        else:
            if validate_subkey_a_size(subkey_a, num_total_qubits):
                subkey_a_mapping(subkey_a,subkey_a_gates,num_total_qubits)
            else:
                validate_subkey_a_size(subkey_a, num_total_qubits, True)
                raise Exception("Invalid key size.")
        #reverse the key so that when the array is passed in as a value to this
        #function and also printed out, it is displayed in reverse order
        #following the convention of Qiskit
        subkey_a_gates.reverse()
    
    subkey_b_permutation = []
    if subkey_b == None:
        m_over_d = int((num_total_qubits - num_signature_qubits)/num_signature_qubits)
        subkey_b_permutation = num_signature_qubits * [num_signature_qubits]
        for i in range(num_signature_qubits):
            subkey_b_permutation.extend(m_over_d * [i])
        random.shuffle(subkey_b_permutation)
    else:
        if validate_subkey_b_size(subkey_b,num_total_qubits, num_signature_qubits):
            subkey_b_permutation = subkey_b_mapping(num_total_qubits, num_signature_qubits, subkey_b)
        else:
            validate_subkey_b_size(subkey_b,num_total_qubits,num_signature_qubits, True)
            raise Exception("Invalid key size.")
            
    #Here, we will say that signature qubit positions are determined by the
    #elements with value 'num_signature_qubits'
    subkey_b_permutation_no_signature_qubits = [x for x in subkey_b_permutation if x != num_signature_qubits]
    
        
    #initialize a quantum circuit for encrypting data
    qc = QuantumCircuit(QuantumRegister(num_total_qubits))
    #apply a sequence of CXs, CYs and CZs, with the controls being the data qubits and the targets being the signature qubits
    for j in range(num_signature_qubits):
        control_gate_to_use = subkey_a_gates[j]
        if control_gate_to_use == 0:
            qc.i(j)
        elif control_gate_to_use == 1:
            qc.h(j)
            my_sdg(qc,j)
        else:
            qc.h(j)
        
        for i in range(num_signature_qubits, num_total_qubits):
            if subkey_b_permutation_no_signature_qubits[i - num_signature_qubits] == j:
                if control_gate_to_use == 0:
                    my_cx(qc,i,j)
                elif control_gate_to_use == 1:
                    my_cz(qc,i,j)
                else:
                    my_cy(qc,i,j)
    
    #apply the unitaries to encrypt the data qubits, determined by some classical encryption key
    #for the purpose of this simulation, the classical encryption key is abstracted away via random unitary generation
    for i in range(num_signature_qubits,num_total_qubits):
        gate_to_use = subkey_a_gates[i]#random int, inclusive for both numbers
        if gate_to_use == 0:
            qc.i(i)
        elif gate_to_use == 1:
            qc.h(i)
            my_sdg(qc,i)
        else:
            qc.h(i)
    
    #move the signature qubits to new positions in the packet based off of
    #subkey_b
    current_j = num_signature_qubits - 1
    for i in range(num_total_qubits-1,-1,-1):
        if subkey_b_permutation[i] == num_signature_qubits:
            my_circ(qc,current_j,i)
            current_j = current_j - 1

    return Clifford(qc)

def block_clifford_decryption(num_total_qubits, num_signature_qubits,subkey_a=None,subkey_b=None,subkey_a_base=2):
    subkey_a_gates = []
    if subkey_a == None:
        for i in range(num_total_qubits):
            subkey_a_gates.append(random.randint(0, 2))
    else:
        if subkey_a_base == 3:
            subkey_a_gates = subkey_a
        else:
            if validate_subkey_a_size(subkey_a, num_total_qubits):
                subkey_a_mapping(subkey_a,subkey_a_gates,num_total_qubits)
            else:
                validate_subkey_a_size(subkey_a, num_total_qubits, True)
                raise Exception("Invalid key size.")
        #reverse the key so that when the array is passed in as a value to this
        #function and also printed out, it is displayed in reverse order
        #following the convention of Qiskit
        subkey_a_gates.reverse()
    
    subkey_b_permutation = []
    if subkey_b == None:
        m_over_d = int((num_total_qubits - num_signature_qubits)/num_signature_qubits)
        subkey_b_permutation = num_signature_qubits * [num_signature_qubits]
        for i in range(num_signature_qubits):
            subkey_b_permutation.extend(m_over_d * [i])
        random.shuffle(subkey_b_permutation)
    else:
        if validate_subkey_b_size(subkey_b,num_total_qubits, num_signature_qubits):
            subkey_b_permutation = subkey_b_mapping(num_total_qubits, num_signature_qubits, subkey_b)
        else:
            validate_subkey_b_size(subkey_b,num_total_qubits,num_signature_qubits, True)
            raise Exception("Invalid key size.")
            
    #Here, we will say that signature qubit positions are determined by the
    #elements with value 'num_signature_qubits'
    subkey_b_permutation_no_signature_qubits = [x for x in subkey_b_permutation if x != num_signature_qubits]
        
    #initialize a quantum circuit for encrypting data
    qc = QuantumCircuit(QuantumRegister(num_total_qubits))
            
    #move the signature qubits to new positions in the packet based off of
    #subkey_b
    current_j = 0
    for i in range(num_total_qubits):
        if subkey_b_permutation[i] == num_signature_qubits:
            my_circ(qc,i,current_j)
            current_j = current_j + 1
    
    #apply the unitaries to encrypt the data qubits, determined by some classical encryption key
    #for the purpose of this simulation, the classical encryption key is abstracted away via random unitary generation
    for i in range(num_signature_qubits,num_total_qubits):
        gate_to_use = subkey_a_gates[i]#random int, inclusive for both numbers
        if gate_to_use == 0:
            qc.i(i)
        elif gate_to_use == 1:
            qc.s(i)
            qc.h(i)
        else:
            qc.h(i)
        
    #apply a sequence of CXs, CYs and CZs, with the controls being the data qubits and the targets being the signature qubits
    for j in range(num_signature_qubits):
        control_gate_to_use = subkey_a_gates[j]
        for i in range(num_signature_qubits, num_total_qubits):
            if subkey_b_permutation_no_signature_qubits[i - num_signature_qubits] == j:
                if control_gate_to_use == 0:
                    my_cx(qc,i,j)
                elif control_gate_to_use == 1:
                    my_cz(qc,i,j)
                else:
                    my_cy(qc,i,j)
                    
        if control_gate_to_use == 0:
            qc.i(j)
        elif control_gate_to_use == 1:
            qc.s(j)
            qc.h(j)
        else:
            qc.h(j)
    

    return Clifford(qc)

def apply_control_gate_from_key_clifford(qc,control,target,control_gate_to_use):
    if control_gate_to_use == 0:
        my_cx(qc,control,target)
    elif control_gate_to_use == 1:
        my_cz(qc,control,target)
    else:
        my_cy(qc,control,target)
        
def apply_single_qubit_gate_from_key_clifford(qc,qubit_position,gate_to_use):
    if gate_to_use == 0:
        qc.i(qubit_position)
    elif gate_to_use == 1:
        qc.h(qubit_position)
        my_sdg(qc,qubit_position)
    else:
        qc.h(qubit_position)

#returns a random state vector for the data qubits, and initializes all signature qubits to |0>
def generate_data_to_encrypt(num_total_qubits, num_signature_bits, interleaved = False):
    signature_qubits = Statevector(np.concatenate((np.ones(1),np.zeros((2 ** num_signature_bits)-1)),axis=0))
    data_qubits = random_statevector(2**(num_total_qubits - num_signature_bits))
    data_to_encrypt = Statevector(np.kron(data_qubits, signature_qubits))
    if interleaved:
        print('interleaved authenticator bits not implemented yet')
    return data_to_encrypt

def adversary_simulation_block_clifford(num_total_qubits, num_signature_qubits, adversarial_attack, data_to_encrypt, subkey_a=None,subkey_b=None,subkey_a_base=2):
    qc = block_clifford_encryption(num_total_qubits, num_signature_qubits,subkey_a=None,subkey_b=None,subkey_a_base=2)
    return adversary_simulation_common(qc, num_total_qubits, num_signature_qubits, adversarial_attack, data_to_encrypt)
    
def adversary_simulation_clifford(num_total_qubits, num_signature_bits, adversarial_attack, data_to_encrypt):
    qc = quantum_info.random_clifford(num_total_qubits)
    return adversary_simulation_common(qc, num_total_qubits, num_signature_bits, adversarial_attack, data_to_encrypt)
    
def adversary_simulation_common(qc, num_total_qubits, num_signature_bits, adversarial_attack, data_to_encrypt, return_transformation=False):    
    if type(qc) is quantum_info.operators.symplectic.clifford.Clifford and type(adversarial_attack) is quantum_info.operators.symplectic.pauli.Pauli:
        resulting_pauli = adversarial_attack.evolve(qc)
        if return_transformation:
            return resulting_pauli
        paulis_affecting_signature_qubits = resulting_pauli.to_label()[-num_signature_bits:]
        if "X" in paulis_affecting_signature_qubits or "Y" in paulis_affecting_signature_qubits:
            return 1
    else:
        
        encrypted_data = data_to_encrypt.evolve(qc)
        attacked_data = encrypted_data.evolve(adversarial_attack)
        decrypted_data = attacked_data.evolve(qc.inverse())
        
            
        measurement_result, resulting_state_vector = decrypted_data.measure(list(range(num_signature_bits)))
        if '1' in str(measurement_result):
            return 1

    return 0