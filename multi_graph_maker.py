import numpy as np
import matplotlib.pyplot as plt
import graph_maker
import utilities
import os
    
def plot_4_dot_graphs():
    fig, axs = plt.subplots(2,2,figsize=(10,10))
    fig.suptitle("Probability of Detecting Different Types of Pauli Attacks")
    ((ax1, ax2), (ax3, ax4)) = axs
    
    main_folder = os.getcwd() + "/Results/Varying_Random_Pauli_Attack/"
    num_experiments = 5000
    
    error_bar_width = 1
    error_bar_cap_size = 3
    dot_size = 3
    
    num_total_qubits = 8
    num_signature_qubits = 4
    num_detections = np.load(main_folder + "m_4_d_4_e_5000_num_detections.npy")
    probability_of_detection = num_detections/num_experiments
    x = list(range(1,probability_of_detection.shape[0]+1))
    
    y = graph_maker.estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)
    ax1.errorbar(x, y, label="Block Clifford code (calculated estimate)",  fmt = "o", markersize=dot_size, color = 'green')

    lower_bound, upper_bound = utilities.wilson_score_confidence(num_detections, num_experiments)
    wilson_score = np.vstack(((lower_bound - probability_of_detection) * -1, upper_bound - probability_of_detection))
    ax1.errorbar(x, probability_of_detection, yerr = wilson_score, fmt = "o", elinewidth = error_bar_width, ecolor = "black", capsize = error_bar_cap_size, capthick = error_bar_width, label="Block Clifford code (simulation)", markersize=dot_size, color = 'blue')
    
    y = graph_maker.calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    ax1.errorbar(x, y, label="Clifford code",  fmt = "o", markersize=dot_size, color = (1,0.4,0))
    
    ax1.grid()
    ax1.set(ylabel='Probability of detection')
    ax1.set_title("n = " + str(num_total_qubits) + ", d = " + str(num_signature_qubits))
    ax1.set_ylim(0.49,1.01)
    ax1.set_xlim(0.5,num_total_qubits + 0.5)
    
    num_total_qubits = 16
    num_signature_qubits = 8
    num_detections = np.load(main_folder + "m_8_d_8_e_5000_num_detections.npy")
    probability_of_detection = num_detections/num_experiments
    x = list(range(1,probability_of_detection.shape[0]+1))

    y = graph_maker.estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)
    ax2.errorbar(x, y, label="Block Clifford code (calculated estimate)",  fmt = "o", markersize=dot_size, color = 'green')
    
    lower_bound, upper_bound = utilities.wilson_score_confidence(num_detections, num_experiments)
    wilson_score = np.vstack(((lower_bound - probability_of_detection) * -1, upper_bound - probability_of_detection))
    ax2.errorbar(x, probability_of_detection, yerr = wilson_score, fmt = "o", elinewidth = error_bar_width, ecolor = "black", capsize = error_bar_cap_size, capthick = error_bar_width, label="Block Clifford code (simulation)", markersize=dot_size, color = 'blue')
    
    y = graph_maker.calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    ax2.errorbar(x, y, label="Clifford code",  fmt = "o", markersize=dot_size, color = (1,0.4,0))
    
    ax2.grid()
    ax2.set_title("n = " + str(num_total_qubits) + ", d = " + str(num_signature_qubits))
    ax2.set_ylim(0.49,1.01)
    ax2.set_xlim(0.5,num_total_qubits + 0.5)
    
    num_total_qubits = 32
    num_signature_qubits = 8
    num_detections = np.load(main_folder + "m_24_d_8_e_5000_num_detections.npy")
    probability_of_detection = num_detections/num_experiments
    x = list(range(1,probability_of_detection.shape[0]+1))
    
    y = graph_maker.estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)
    ax3.errorbar(x, y, label="Block Clifford code (calculated estimate)",  fmt = "o", markersize=dot_size, color = 'green')
    
    lower_bound, upper_bound = utilities.wilson_score_confidence(num_detections, num_experiments)
    wilson_score = np.vstack(((lower_bound - probability_of_detection) * -1, upper_bound - probability_of_detection))
    ax3.errorbar(x, probability_of_detection, yerr = wilson_score, fmt = "o", elinewidth = error_bar_width, ecolor = "black", capsize = error_bar_cap_size, capthick = error_bar_width, label="Block Clifford code (simulation)", markersize=dot_size, color = 'blue')
    
    y = graph_maker.calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    ax3.errorbar(x, y, label="Clifford code",  fmt = "o", markersize=dot_size, color = (1,0.4,0))
    
    ax3.grid()
    ax3.set(xlabel='Number of non-identity single qubit Paulis')
    ax3.set(ylabel='Probability of detection')
    ax3.set_title("n = " + str(num_total_qubits) + ", d = " + str(num_signature_qubits))
    ax3.set_ylim(0.49,1.01)
    ax3.set_xlim(0.5,num_total_qubits + 0.5)
    
    num_total_qubits = 32
    num_signature_qubits = 16
    num_detections = np.load(main_folder + "m_16_d_16_e_5000_num_detections.npy")

    probability_of_detection = num_detections/num_experiments
    x = list(range(1,probability_of_detection.shape[0]+1))
    
    y = graph_maker.estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)
    ax4.errorbar(x, y, label="Block Clifford code (calculated estimate)",  fmt = "o", markersize=dot_size, color = 'green')

    lower_bound, upper_bound = utilities.wilson_score_confidence(num_detections, num_experiments)
    wilson_score = np.vstack(((lower_bound - probability_of_detection) * -1, upper_bound - probability_of_detection))
    ax4.errorbar(x, probability_of_detection, yerr = wilson_score, fmt = "o", elinewidth = error_bar_width, ecolor = "black", capsize = error_bar_cap_size, capthick = error_bar_width, label="Block Clifford code (simulation)", markersize=dot_size, color = 'blue')
    
    y = graph_maker.calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    ax4.errorbar(x, y, label="Clifford code",  fmt = "o", markersize=dot_size, color = (1,0.4,0))
    
    ax4.set(xlabel='Number of non-identity single qubit Paulis')
    ax4.grid()
    ax4.legend()
    ax4.set_title("n = " + str(num_total_qubits) + ", d = " + str(num_signature_qubits))
    ax4.set_ylim(0.49,1.01)
    ax4.set_xlim(0.5,num_total_qubits + 0.5)
    
    fig.tight_layout()
    plt.savefig('multi_dot_graph_pauli_attacks.png', format='png', dpi=600)
    plt.show()

# plot_4_dot_graphs()

def plot_4_line_graphs():
    fig, axs = plt.subplots(2,2,figsize=(10,10))
    fig.suptitle("Probability of Detecting Different Types of Pauli Attacks")
    ((ax1, ax2), (ax3, ax4)) = axs
    
    main_folder = os.getcwd() + "/Results/Varying_Random_Pauli_Attack/"
    num_experiments = 5000
    
    num_total_qubits = 8
    num_signature_qubits = 4
    num_detections = np.load(main_folder + "m_4_d_4_e_5000_num_detections.npy")
    probability_of_detection = num_detections/num_experiments
    x = list(range(1,probability_of_detection.shape[0]+1))
    lower_bound, upper_bound = utilities.wilson_score_confidence(num_detections, num_experiments)
    ax1.plot(x, probability_of_detection, label="Block Clifford code (simulation)", linewidth=0.65)
    ax1.fill_between(x, lower_bound, upper_bound, alpha=.1)
    y = graph_maker.calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    ax1.plot(x, y, label="Clifford code",  linewidth=0.65)
    y = graph_maker.estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)
    ax1.plot(x, y, label="Block Clifford code (calculated estimate)",  linewidth=0.65)
    ax1.set(ylabel='Probability of detection')
    ax1.set_title("n = " + str(num_total_qubits) + ", d = " + str(num_signature_qubits))
    ax1.set_ylim(0.49,1.01)
    ax1.set_xlim(1,num_total_qubits)
    
    num_total_qubits = 16
    num_signature_qubits = 8
    num_detections = np.load(main_folder + "m_8_d_8_e_5000_num_detections.npy")
    probability_of_detection = num_detections/num_experiments
    x = list(range(1,probability_of_detection.shape[0]+1))
    lower_bound, upper_bound = utilities.wilson_score_confidence(num_detections, num_experiments)
    ax2.plot(x, probability_of_detection, label="Block Clifford code (simulation)", linewidth=0.65)
    ax2.fill_between(x, lower_bound, upper_bound, alpha=.1)
    y = graph_maker.calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    ax2.plot(x, y, label="Clifford code",  linewidth=0.65)
    y = graph_maker.estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)
    ax2.plot(x, y, label="Block Clifford code (calculated estimate)",  linewidth=0.65)
    ax2.set_title("n = " + str(num_total_qubits) + ", d = " + str(num_signature_qubits))
    ax2.set_ylim(0.49,1.01)
    ax2.set_xlim(1,num_total_qubits)
    
    num_total_qubits = 32
    num_signature_qubits = 8
    num_detections = np.load(main_folder + "m_24_d_8_e_5000_num_detections.npy")
    probability_of_detection = num_detections/num_experiments
    x = list(range(1,probability_of_detection.shape[0]+1))
    lower_bound, upper_bound = utilities.wilson_score_confidence(num_detections, num_experiments)
    ax3.plot(x, probability_of_detection, label="Block Clifford code (simulation)", linewidth=0.65)
    ax3.fill_between(x, lower_bound, upper_bound, alpha=.1)
    y = graph_maker.calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    ax3.plot(x, y, label="Clifford code",  linewidth=0.65)
    y = graph_maker.estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)
    ax3.plot(x, y, label="Block Clifford code (calculated estimate)",  linewidth=0.65)
    ax3.set(xlabel='Number of non-identity single qubit Paulis', ylabel='Probability of detection')
    ax3.set_title("n = " + str(num_total_qubits) + ", d = " + str(num_signature_qubits))
    ax3.set_ylim(0.49,1.01)
    ax3.set_xlim(1,num_total_qubits)
    
    num_total_qubits = 32
    num_signature_qubits = 16
    num_detections = np.load(main_folder + "m_16_d_16_e_5000_num_detections.npy")

    probability_of_detection = num_detections/num_experiments
    x = list(range(1,probability_of_detection.shape[0]+1))
    lower_bound, upper_bound = utilities.wilson_score_confidence(num_detections, num_experiments)
    ax4.plot(x, probability_of_detection, label="Block Clifford code (simulation)", linewidth=0.65)
    ax4.fill_between(x, lower_bound, upper_bound, alpha=.1)
    y = graph_maker.calculate_clifford_code_detection_rate(num_total_qubits,num_signature_qubits)
    ax4.plot(x, y, label="Clifford code",  linewidth=0.65)
    y = graph_maker.estimate_block_clifford_code_detection_rate(num_signature_qubits, num_total_qubits)

    ax4.plot(x, y, label="Block Clifford code (calculated estimate)",  linewidth=0.65)
    ax4.set(xlabel='Number of non-identity single qubit Paulis')
    ax4.legend()
    ax4.set_title("n = " + str(num_total_qubits) + ", d = " + str(num_signature_qubits))
    ax4.set_ylim(0.49,1.01)
    ax4.set_xlim(1,num_total_qubits)
    
    fig.tight_layout()
    plt.savefig('multi_line_graph_pauli_attacks.png', format='png', dpi=600)
    plt.show()
    
# plot_4_line_graphs()