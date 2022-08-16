import numpy as np
import math
from qiskit import QuantumRegister
from qiskit import QuantumCircuit
from qiskit.tools.visualization import circuit_drawer
import encryption


#https://github.com/NickPerezCarletonUniversity/Authenticity-Integrity-and-Replay-Protection-in-Quantum-Data-Communications-and-Networking/blob/main/utilities.py
#https://stackoverflow.com/questions/10029588/python-implementation-of-the-wilson-score-interval
#https://www.dummies.com/education/math/statistics/checking-out-statistical-confidence-interval-critical-values/
def wilson_score_confidence(successes, trials, z=1.96):
    z = 3.08
    n = trials

    if n == 0:
        return 0

    upper_bound = np.zeros(successes.shape[0])
    lower_bound = np.zeros(successes.shape[0])
    
    for i in range(successes.shape[0]):
        phat = float(successes[i]) / n
        upper_bound[i] = ((phat + z*z/(2*n) + z * math.sqrt((phat*(1-phat)+z*z/(4*n))/n))/(1+z*z/n))
        lower_bound[i] = ((phat + z*z/(2*n) - z * math.sqrt((phat*(1-phat)+z*z/(4*n))/n))/(1+z*z/n))
    return lower_bound, upper_bound


def create_example_circuit_diagram():
    q = QuantumRegister(4)
    qc = QuantumCircuit(q)
    qc.h(0)
    qc.cy(2,0)
    qc.barrier()
    qc.i(1)
    qc.cx(3,1)
    qc.barrier()
    qc.h(2)
    qc.barrier()
    qc.h(3)
    qc.sdg(3)
    qc.barrier()
    encryption.my_circ_simple(qc, 0, 2)
    
    print(qc)
    circuit_drawer(qc, output='mpl',filename="test_circuit", style={'backgroundcolor': '#EEEEEE'})


# create_example_circuit_diagram()