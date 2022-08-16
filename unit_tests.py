import encryption
import random
from qiskit.quantum_info import Pauli
from qiskit import QuantumRegister
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator
import math

def determine_max_key_size(num_total_qubits,num_signature_qubits):
    return (encryption.determine_max_subkey_a(num_total_qubits)[1] + 1) * (encryption.determine_max_subkey_b(num_total_qubits, num_signature_qubits)[1] + 1)

# print(len(str(determine_max_key_size(32,16))))
# print(determine_max_key_size(32,16))

def verify_subkey_b_mapping():
    for d_multiplier in range(3,4):
        for d in range(3,4):
            n = d * d_multiplier
            combinations = set()
            num_combinations = 1
            m_over_d = int((n - d)/d)
            for i in range(d):
                num_combinations *= math.comb(n - i * m_over_d, m_over_d)
            for subkey in range(num_combinations):
                combination = encryption.subkey_b_mapping(n,d,subkey,False)
                combination_str = str(combination)
                overflow_combination = encryption.subkey_b_mapping(n,d,subkey + num_combinations,False)
                overflow_combination_str = str(overflow_combination)
                if overflow_combination_str != combination_str:
                    return False
                for i in range(d):
                    if combination.count(i) != m_over_d:
                        return False
                if combination_str in combinations or len(combination) != n or combination.count(d) != d:
                    return False
                else:
                    combinations.add(combination_str)
    return True

print(verify_subkey_b_mapping())

def verify_quick_n_choose_k():
    prev_n = None
    prev_k = None
    n_choose_k = 0
    for k in range(1,5):
        n_choose_k = 0
        for n in range(k,10):
            n_choose_k, prev_n, prev_k = encryption.quick_n_choose_k(n, k, prev_n, prev_k, n_choose_k)
            if n_choose_k != math.comb(n,k):
                return False
    for k in range(1,5):
        n_choose_k = math.comb(11,k)
        prev_n = 11
        prev_k = k
        for n in range(10,k - 1,-1):
            n_choose_k, prev_n, prev_k = encryption.quick_n_choose_k(n, k, prev_n, prev_k, n_choose_k)
            if n_choose_k != math.comb(n,k):
                return False
    for k in range(0,5):
        n_choose_k = math.comb(11,k)
        prev_n = 11
        prev_k = 11-k
        for n in range(10,5,-1):
            n_choose_k, prev_n, prev_k = encryption.quick_n_choose_k(n, n-k, prev_n, prev_k, n_choose_k)
            if n_choose_k != math.comb(n,k):
                return False
    return True

print(verify_quick_n_choose_k())

def verify_combinatorial_number_system():
    for n in range(12):
        for k in range(1,n+1):
            combinations = set()
            for N in range(math.comb(n,k)):
                combination = encryption.combinatorial_number_system(n,k,N)
                combination_str = str(combination)
                for element in combination:
                    if element < 0 or element >= n:
                        return False
                if combination_str in combinations or len(combination) != k or len(set(combination)) != len(combination):
                    return False
                else:
                    combinations.add(combination_str)
    return True

print(verify_combinatorial_number_system())

def verify_block_clifford_code_inverse():
    num_total_qubits = 8
    num_signature_qubits = 2
    subkey_a = []
    subkey_b = []
    for i in range(encryption.determine_max_subkey_a(num_total_qubits)[0]):
        subkey_a.append(random.randint(0,1))
    for i in range(encryption.determine_max_subkey_b(num_total_qubits,num_signature_qubits)[0]):
        subkey_b.append(random.randint(0,1))
    qc_encrypt = encryption.block_clifford_encryption(num_total_qubits, num_signature_qubits,subkey_a,subkey_b)
    qc_decrypt = encryption.block_clifford_decryption(num_total_qubits, num_signature_qubits,subkey_a,subkey_b)
    if not Operator(qc_encrypt.adjoint()).equiv(Operator(qc_decrypt)):
        return False
    return True

print(verify_block_clifford_code_inverse())


#verifies detection probabilities for n = 2 and d = 1
def verify_two_qubit_pauli_attack_theorem_with_only_two_qubits():
    num_total_qubits = 2
    num_signature_qubits = 1
    num_experiments = 2000
     
    adversarial_attacks = [Pauli("IX"), Pauli("IY"), Pauli("IZ"), Pauli("XX"), Pauli("XY"), Pauli("XZ"), 
                           Pauli("YX"), Pauli("YY"), Pauli("YZ"), Pauli("ZX"), Pauli("ZY"), Pauli("ZZ"),
                           Pauli("XI"), Pauli("YI"), Pauli("ZI")]

    for pauli_attack in adversarial_attacks:
        num_detections = 0
        for experiment in range(num_experiments):
            if encryption.adversary_simulation_block_clifford(num_total_qubits, num_signature_qubits, pauli_attack, None) == 1:
                num_detections = num_detections + 1
        print(num_detections/num_experiments)
        
    return True

# print(verify_two_qubit_pauli_attack_theorem_with_only_two_qubits())

#shows all the different resulting transformations from a pauli attack on a
#data qubit with only one data qubit and one signature qubit in the packet
#subkey B is set to the identity permutation
def verify_single_qubit_pauli_attack_on_a_signature_qubit_theorem_with_only_two_qubits(verbose = False):
    num_total_qubits = 2
    num_signature_qubits = 1
     
    adversarial_attacks = [Pauli("IX"), Pauli("IY"), Pauli("IZ")]
    
    for signature_qubit_encrypting_gate in range(3):
        for data_qubit_encrypting_gate in range(3):
            for pauli_attack in adversarial_attacks:
                subkey_a = [data_qubit_encrypting_gate, signature_qubit_encrypting_gate]
                qc = encryption.block_clifford_encryption(num_total_qubits, num_signature_qubits, subkey_a, [1], 3)
                resulting_pauli = pauli_attack.evolve(qc)
                if verbose:
                    print("When using the subkey A (in base 3): " + str([data_qubit_encrypting_gate, signature_qubit_encrypting_gate]))
                    print("If the adversary attacks using:      " + str(pauli_attack))
                    print("The resulting Pauli is:              " + str(resulting_pauli) + "\n")

    return True

# print(verify_single_qubit_pauli_attack_on_a_signature_qubit_theorem_with_only_two_qubits(True))

#shows all the different resulting transformations from a pauli attack on a
#data qubit with only one data qubit and one signature qubit in the packet
#subkey B is set to the identity permutation
def verify_single_qubit_pauli_attack_on_a_data_qubit_theorem_with_only_two_qubits(verbose = False):
    num_total_qubits = 2
    num_signature_qubits = 1
     
    adversarial_attacks = [Pauli("XI"), Pauli("YI"), Pauli("ZI")]
    
    for signature_qubit_encrypting_gate in range(3):
        for data_qubit_encrypting_gate in range(3):
            for pauli_attack in adversarial_attacks:
                subkey_a = [data_qubit_encrypting_gate, signature_qubit_encrypting_gate]
                qc = encryption.block_clifford_encryption(num_total_qubits, num_signature_qubits, subkey_a, [1], 3)
                resulting_pauli = pauli_attack.evolve(qc)
                if verbose:
                    print("When using the subkey A (in base 3): " + str([data_qubit_encrypting_gate, signature_qubit_encrypting_gate]))
                    print("If the adversary attacks using:      " + str(pauli_attack))
                    print("The resulting Pauli is:              " + str(resulting_pauli) + "\n")

    return True

# print(verify_single_qubit_pauli_attack_on_a_data_qubit_theorem_with_only_two_qubits(True))


def verify_my_sdg():
    num_qubits = 6
    for target in range(num_qubits):
        qc1 = QuantumCircuit(QuantumRegister(num_qubits))
        qc1.sdg(target)

        qc2 = QuantumCircuit(QuantumRegister(num_qubits))
        encryption.my_sdg(qc2,target)

        if not Operator(qc1).equiv(Operator(qc2)):
            return False
    return True
    
print(verify_my_sdg())

def verify_my_swap():
    num_qubits = 6
    for target in range(num_qubits-1):
        qc1 = QuantumCircuit(QuantumRegister(num_qubits))
        qc1.swap(target,target+1)

        qc2 = QuantumCircuit(QuantumRegister(num_qubits))
        encryption.my_swap(qc2,target)

        if not Operator(qc1).equiv(Operator(qc2)):
            return False
    return True
    
print(verify_my_swap())

def verify_my_circ():
    num_qubits = 6
    for source in range(num_qubits):
        for destination in range(num_qubits):
            if source != destination:
                qc1 = QuantumCircuit(QuantumRegister(num_qubits))
                encryption.my_circ(qc1,source,destination)

                qc2 = QuantumCircuit(QuantumRegister(num_qubits))
                encryption.my_circ(qc2,destination,source)

                if not Operator(qc1.inverse()).equiv(Operator(qc2)):
                    return False
    return True
            
print(verify_my_circ())

def verify_my_cx():
    num_qubits = 6
    for control in range(num_qubits):
        for target in range(num_qubits):
            if control != target:
                qc1 = QuantumCircuit(QuantumRegister(num_qubits))
                qc1.cx(control,target)

                qc2 = QuantumCircuit(QuantumRegister(num_qubits))
                encryption.my_cx(qc2,control,target)

                if not Operator(qc1).equiv(Operator(qc2)):
                    return False
    return True
            
print(verify_my_cx())

def verify_my_cy():
    num_qubits = 6
    for control in range(num_qubits):
        for target in range(num_qubits):
            if control != target:
                qc1 = QuantumCircuit(QuantumRegister(num_qubits))
                qc1.cy(control,target)
                
                qc2 = QuantumCircuit(QuantumRegister(num_qubits))
                encryption.my_cy(qc2,control,target)

                if not Operator(qc1).equiv(Operator(qc2)):
                    return False
    return True
            
print(verify_my_cy())

def verify_my_cz():
    num_qubits = 6
    for control in range(num_qubits):
        for target in range(num_qubits):
            if control != target:
                qc1 = QuantumCircuit(QuantumRegister(num_qubits))
                qc1.cz(control,target)

                qc2 = QuantumCircuit(QuantumRegister(num_qubits))
                encryption.my_cz(qc2,control,target)

                if not Operator(qc1).equiv(Operator(qc2)):
                    return False
    return True
            
print(verify_my_cz())

def verify_subkey_a_mapping():
    num_qubits = 4
    key = [1,0,0,1,1,0,1]
    gates_to_use = []
    encryption.subkey_a_mapping(key,gates_to_use,num_qubits)
    if gates_to_use != [2, 2, 1, 2]:
        return False

    key = [0,1,1,0,1,1]
    gates_to_use = []
    encryption.subkey_a_mapping(key,gates_to_use,num_qubits)
    if gates_to_use != [1,0,0,0]:
        return False
    
    key = [0,0,0,0,0,1]
    gates_to_use = []
    encryption.subkey_a_mapping(key,gates_to_use,num_qubits)
    if gates_to_use != [0,0,0,1]:
        return False
    
    return True
    
print(verify_subkey_a_mapping())

def verify_validate_subkey_a_size():
    key = [0,0,0,0,0,0,1]
    num_qubits = 4
    if not encryption.validate_subkey_a_size(key,num_qubits):
        return False
    
    num_qubits = 10
    if encryption.validate_subkey_a_size(key,num_qubits):
        return False
    
    key = [0,0,0,0,0,0,0,0]
    num_qubits = 5
    if not encryption.validate_subkey_a_size(key,num_qubits):
        return False
    
    num_qubits = 1
    if encryption.validate_subkey_a_size(key,num_qubits):
        return False
    
    key = [1,0,0,1,1,0,0,1,1,1,0,0,0,1,0]
    num_qubits = 9
    if not encryption.validate_subkey_a_size(key,num_qubits):
        return False
    
    num_qubits = 11
    if encryption.validate_subkey_a_size(key,num_qubits):
        return False
    
    return True
    
print(verify_validate_subkey_a_size())

