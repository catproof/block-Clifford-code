import numpy as np
import time
import adversarial_attack as attacks
import encryption
import matplotlib.pyplot as plt
import math
import utilities
import os

def calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits):
    clifford_code_detection = np.ones(num_total_qubits)
    clifford_code_detection_equation = 1 - (2**(2 * (num_total_qubits - num_signature_qubits) + num_signature_qubits) - 1)/(2**(2 * num_total_qubits) - 1)
    clifford_code_detection = clifford_code_detection * clifford_code_detection_equation
    return clifford_code_detection

def estimate_block_clifford_code_detection_rate_helper(num_attacks, num_signature_qubits):
    w = num_attacks
    d = num_signature_qubits
    probability_of_detection = 0
    for x in range(1,min(w,d) + 1):
        probability_of_num_attacked_blocks = 0
        for i in range(x+1):
            probability_of_num_attacked_blocks += (-1)**i * math.comb(x,i) * ((x-i)/d)**w
        probability_of_num_attacked_blocks *= math.comb(d,x)
        probability_of_detection += (1 - (1/2)**x) * probability_of_num_attacked_blocks
    return probability_of_detection

def verify_estimate_block_clifford_code_detection_rate_helper():
    for num_attacks in range(1,64):
        for num_signature_qubits in range(1,32):
            w = num_attacks
            d = num_signature_qubits
            sum_of_probabilities = 0
            for x in range(1,min(w,d) + 1):
                probability_of_num_attacked_blocks = 0
                for i in range(x+1):
                    probability_of_num_attacked_blocks += (-1)**i * math.comb(x,i) * ((x-i)/d)**w
                probability_of_num_attacked_blocks *= math.comb(d,x)
                sum_of_probabilities = sum_of_probabilities + probability_of_num_attacked_blocks
            if abs(sum_of_probabilities - 1) > 0.0000001:
                print(sum_of_probabilities)
                return False
            
    return True

# print(verify_estimate_block_clifford_code_detection_rate_helper())

def estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits):
    block_clifford_code_detection_rate_estimate = np.zeros(num_total_qubits)
    for num_non_identity_paulis in range(1, num_total_qubits + 1):
        block_clifford_code_detection_rate_estimate[num_non_identity_paulis - 1] = estimate_block_clifford_code_detection_rate_helper(num_non_identity_paulis, num_signature_qubits)
    return block_clifford_code_detection_rate_estimate

# print(estimate_block_clifford_code_detection_rate_helper(20, 8))
# num_total_qubits = 64
# num_signature_qubits = 16
# print(estimate_block_clifford_code_detection_rate(num_total_qubits, num_signature_qubits))
# print(calculate_clifford_code_detection_rate(num_total_qubits, num_signature_qubits))

def block_clifford_code_experiment_paulis(num_trials,num_total_qubits,num_signature_qubits):
    
    start = time.time()
    
    adversarial_attack = attacks.specific_random_paulis(num_total_qubits,min_non_identity_paulis = num_total_qubits, max_non_identity_paulis = num_total_qubits)
    encryption.adversary_simulation_block_clifford(num_total_qubits, num_signature_qubits, adversarial_attack, None)
    
    stop = time.time()
    print("estimated execution time of experiment: " + str((stop - start) * num_trials * num_total_qubits / 60) + " minutes.")
    
    print("n = " + str(num_total_qubits))
    print("d = " + str(num_signature_qubits))
    
    block_clifford_code_num_detections_list = np.zeros(num_total_qubits)
    for num_non_identity_paulis in range(1, num_total_qubits + 1):
        print("Running trials for Pauli attacks with " + str(num_non_identity_paulis) + " non-identity single-qubit Paulis")
        for experiment in range(num_trials):
    
            #The adversary applies a Pauli attack
            adversarial_attack = attacks.specific_random_paulis(num_total_qubits,min_non_identity_paulis = num_non_identity_paulis, max_non_identity_paulis = num_non_identity_paulis)
            
            if encryption.adversary_simulation_block_clifford(num_total_qubits, num_signature_qubits, adversarial_attack, None) == 1:
                block_clifford_code_num_detections_list[num_non_identity_paulis - 1] = block_clifford_code_num_detections_list[num_non_identity_paulis - 1] + 1
    
    stop = time.time()
    print("execution time of experiment: " + str((stop - start)/60) + " minutes.")
    
    block_clifford_code_detection = block_clifford_code_num_detections_list/num_trials
    clifford_code_detection = calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    block_clifford_code_detection_calculated_estimate = estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)
    
    block_clifford_code_detection_filename = "m_" + str(num_total_qubits - num_signature_qubits) + "_d_" + str(num_signature_qubits) + "_e_" + str(num_trials)
    np.save(block_clifford_code_detection_filename + str("_detection_rate"),block_clifford_code_detection)
    np.save(block_clifford_code_detection_filename + str("_num_detections"),block_clifford_code_num_detections_list)
    
    return block_clifford_code_detection, block_clifford_code_detection_filename, clifford_code_detection, block_clifford_code_detection_calculated_estimate

# num_trials = 1
# num_total_qubits = 32
# num_signature_qubits = 8

# block_clifford_code_detection, block_clifford_code_detection_filename, clifford_code_detection, block_clifford_code_detection_calculated_estimate = block_clifford_code_experiment_paulis(num_trials,num_total_qubits,num_signature_qubits)

# fig, ax = plt.subplots()
# x = np.arange(num_total_qubits) + 1
# ax.plot(x, block_clifford_code_detection, label='Block Clifford code (simulation)',linewidth=0.5)
# ax.plot(x, clifford_code_detection, label='Clifford code',linewidth=0.5)
# ax.plot(x, block_clifford_code_detection_calculated_estimate, label='Block Clifford code (calculated estimate)',linewidth=0.5)

# ax.set(xlabel='Number of Non-Identity Paulis', ylabel='Probability of Detection', title='Probability of Detecting Different Types of Adversaries')

# ax.legend()
# fig.savefig(block_clifford_code_detection_filename + ".png", format='png', dpi=300)
# plt.show()

def clifford_key_mapping_complexity(n):
    return n**3

def block_clifford_key_mapping_complexity(n):
    return n**2 * (math.log2(n))**2

def key_mapping_time_complexity_comparison(max_num_total_qubits = 4608):
    clifford_key_mapping_time = np.zeros(max_num_total_qubits)
    block_clifford_key_mapping_time = np.zeros(max_num_total_qubits)
    for i in range(1,max_num_total_qubits):
        clifford_key_mapping_time[i] = clifford_key_mapping_complexity(i)
        block_clifford_key_mapping_time[i] = block_clifford_key_mapping_complexity(i)
        
    fig, ax = plt.subplots(figsize=(7,7))
    x = np.arange(max_num_total_qubits) + 1
    ax.plot(x, clifford_key_mapping_time, label='Clifford code',linewidth=0.65, color = (1,0.4,0))
    ax.plot(x, block_clifford_key_mapping_time, label='Block Clifford code',linewidth=0.65, color = 'blue')

    ax.set(xlabel='Number of qubits in the data packet', ylabel='Asymptotic worst-case run time', title='Asymptotic Run Times of Key Mappings')

    ax.legend()
    fig.savefig("keymapping_runtimes.png", format='png', dpi=600)
    plt.show()
    
key_mapping_time_complexity_comparison()

def clifford_gates_complexity(n):
    return n**2 / math.log2(n)

def block_clifford_gates_complexity(n, d):
    return n*d

def gates_complexity_comparison(max_num_total_qubits = 4608, num_signature_qubits = 64):
    clifford_num_gates = np.zeros(max_num_total_qubits)
    block_clifford_num_gates = np.zeros(max_num_total_qubits)
    for i in range(2,max_num_total_qubits):
        clifford_num_gates[i] = clifford_gates_complexity(i)
        block_clifford_num_gates[i] = block_clifford_gates_complexity(i,num_signature_qubits)
        
    fig, ax = plt.subplots(figsize=(7,7))
    x = np.arange(max_num_total_qubits) + 2
    ax.plot(x, clifford_num_gates, label='Clifford code',linewidth=0.65, color = (1,0.4,0))
    ax.plot(x, block_clifford_num_gates, label='Block Clifford code (with d = 64)',linewidth=0.65, color = 'blue')

    ax.set(xlabel='Number of qubits in the data packet', ylabel='Asymptotic worst-case number of gates', title='Asymptotic Number of Gates Used')

    ax.legend()
    fig.savefig("number_of_gates.png", format='png', dpi=600)
    plt.show()
    
# gates_complexity_comparison()

def detection_accuracy_dot_graph(num_total_qubits, num_signature_qubits, file_name, num_experiments = 5000):
    
    fig, ax = plt.subplots(figsize=(7,7))
    
    main_folder = os.getcwd() + "/Results/Varying_Random_Pauli_Attack/"
    num_experiments = 5000
    
    num_detections = np.load(main_folder + file_name)
    probability_of_detection = num_detections/num_experiments
    x = np.arange(1,len(probability_of_detection)+1)
    
    error_bar_width = 1
    error_bar_cap_size = 3
    dot_size = 3
    
    y = estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)
    ax.errorbar(x, y, label="Block Clifford code (calculated estimate)",  fmt = "o", markersize=dot_size, color = 'green')
    
    lower_bound, upper_bound = utilities.wilson_score_confidence(num_detections, num_experiments)
  
    wilson_score = np.vstack(((lower_bound - probability_of_detection) * -1, upper_bound - probability_of_detection))
    
    ax.errorbar(x, probability_of_detection, yerr = wilson_score, fmt = "o", elinewidth = error_bar_width, ecolor = "black", color = 'blue', capsize = error_bar_cap_size, capthick = error_bar_width, label="Block Clifford code (simulation)", markersize=dot_size)
    
    y = calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    ax.errorbar(x, y, label="Clifford code",  fmt = "o", markersize=dot_size, color = (1,0.4,0))
    
    ax.set_ylim(0.49,1.01)
    ax.set_xlim(0.5,num_total_qubits + 0.5)

    ax.set(xlabel='Number of on-identity single-qubit Paulis', ylabel='Probability of detection', title='Probability of Detecting Different Types of Pauli Attacks, n = ' + str(num_total_qubits) + ', d = ' + str(num_signature_qubits))

    ax.grid()    

    ax.legend(loc = 'lower right')
    fig.savefig(file_name[:-4] + "_dot_graph.png", format='png', dpi=600)
    plt.show()
    
# detection_accuracy_dot_graph(8, 4, "m_4_d_4_e_5000_num_detections.npy", num_experiments = 5000)
# detection_accuracy_dot_graph(16, 8, "m_8_d_8_e_5000_num_detections.npy", num_experiments = 5000)
# detection_accuracy_dot_graph(32, 8, "m_24_d_8_e_5000_num_detections.npy", num_experiments = 5000)
# detection_accuracy_dot_graph(32, 16, "m_16_d_16_e_5000_num_detections.npy", num_experiments = 5000)

def detection_accuracy_line_graph(num_total_qubits, num_signature_qubits, file_name, num_experiments = 5000):
    
    fig, ax = plt.subplots(figsize=(7,7))
    
    main_folder = os.getcwd() + "/Results/Varying_Random_Pauli_Attack/"
    num_experiments = 5000
    
    num_detections = np.load(main_folder + file_name)
    probability_of_detection = num_detections/num_experiments
    x = list(range(1,probability_of_detection.shape[0]+1))
    lower_bound, upper_bound = utilities.wilson_score_confidence(num_detections, num_experiments)
    ax.plot(x, probability_of_detection, label="Block Clifford code (simulation)", linewidth=0.65)
    ax.fill_between(x, lower_bound, upper_bound, alpha=.1)
    y = calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    ax.plot(x, y, label="Clifford code",  linewidth=0.65)
    y = estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)
    ax.plot(x, y, label="Block Clifford code (calculated estimate)",  linewidth=0.65)
    ax.set_ylim(0.5,1.01)
    ax.set_xlim(1,num_total_qubits)

    ax.set(xlabel='Number of non-identity single-qubit Paulis', ylabel='Probability of detection', title='Probability of Detecting Different Types of Pauli Attacks, n = ' + str(num_total_qubits) + ', d = ' + str(num_signature_qubits))

    ax.legend(loc = 'lower right')
    fig.savefig(file_name[:-4] + "_line_graph.png", format='png', dpi=600)
    plt.show()

# detection_accuracy_line_graph(8, 4, "m_4_d_4_e_5000_num_detections.npy", num_experiments = 5000)
# detection_accuracy_line_graph(16, 8, "m_8_d_8_e_5000_num_detections.npy", num_experiments = 5000)
# detection_accuracy_line_graph(32, 8, "m_24_d_8_e_5000_num_detections.npy", num_experiments = 5000)
# detection_accuracy_line_graph(32, 16, "m_16_d_16_e_5000_num_detections.npy", num_experiments = 5000)