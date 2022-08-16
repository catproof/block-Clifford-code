from qiskit.quantum_info import Pauli
from qiskit import QuantumRegister
from qiskit import QuantumCircuit


def hadamard_gate():
    qc = QuantumCircuit(QuantumRegister(1))
    qc.h(0)
    qc = qc.inverse()
    print("HXH = " + Pauli("X").evolve(qc).to_label())
    print("HZH = " + Pauli("Z").evolve(qc).to_label())
    print("")
    
hadamard_gate()


def s_gate():
    qc = QuantumCircuit(QuantumRegister(1))
    qc.s(0)
    qc = qc.inverse()
    print("SXS† = " + Pauli("X").evolve(qc).to_label())
    print("SZS† = " + Pauli("Z").evolve(qc).to_label())
    print("")
    
s_gate()

def cx_gate():
    qc = QuantumCircuit(QuantumRegister(2))
    qc.cx(1,0)
    qc = qc.inverse()
    print("CX(X ⊗ I)CX† = " + Pauli("XI").evolve(qc).to_label()[0] + " ⊗ " + Pauli("XI").evolve(qc).to_label()[1])
    print("CX(I ⊗ X)CX† = " + Pauli("IX").evolve(qc).to_label()[0] + " ⊗ " + Pauli("IX").evolve(qc).to_label()[1])
    print("CX(Z ⊗ I)CX† = " + Pauli("ZI").evolve(qc).to_label()[0] + " ⊗ " + Pauli("ZI").evolve(qc).to_label()[1])
    print("CX(I ⊗ Z)CX† = " + Pauli("IZ").evolve(qc).to_label()[0] + " ⊗ " + Pauli("IZ").evolve(qc).to_label()[1])
    print("")
    
cx_gate()

def cy_gate():
    qc = QuantumCircuit(QuantumRegister(2))
    qc.cy(1,0)
    qc = qc.inverse()
    print("CY(X ⊗ I)CY† = " + Pauli("XI").evolve(qc).to_label()[0] + " ⊗ " + Pauli("XI").evolve(qc).to_label()[1])
    print("CY(I ⊗ X)CY† = " + Pauli("IX").evolve(qc).to_label()[0] + " ⊗ " + Pauli("IX").evolve(qc).to_label()[1])
    print("CY(Z ⊗ I)CY† = " + Pauli("ZI").evolve(qc).to_label()[0] + " ⊗ " + Pauli("ZI").evolve(qc).to_label()[1])
    print("CY(I ⊗ Z)CY† = " + Pauli("IZ").evolve(qc).to_label()[0] + " ⊗ " + Pauli("IZ").evolve(qc).to_label()[1])
    print("")
    
cy_gate()

def cz_gate():
    qc = QuantumCircuit(QuantumRegister(2))
    qc.cz(1,0)
    qc = qc.inverse()
    print("CZ(X ⊗ I)CZ† = " + Pauli("XI").evolve(qc).to_label()[0] + " ⊗ " + Pauli("XI").evolve(qc).to_label()[1])
    print("CZ(I ⊗ X)CZ† = " + Pauli("IX").evolve(qc).to_label()[0] + " ⊗ " + Pauli("IX").evolve(qc).to_label()[1])
    print("CZ(Z ⊗ I)CZ† = " + Pauli("ZI").evolve(qc).to_label()[0] + " ⊗ " + Pauli("ZI").evolve(qc).to_label()[1])
    print("CZ(I ⊗ Z)CZ† = " + Pauli("IZ").evolve(qc).to_label()[0] + " ⊗ " + Pauli("IZ").evolve(qc).to_label()[1])
    print("")
    
cz_gate()

def swap_gate():
    qc = QuantumCircuit(QuantumRegister(2))
    qc.swap(1,0)
    qc = qc.inverse()
    print("SWAP(X ⊗ I)SWAP† = " + Pauli("XI").evolve(qc).to_label()[0] + " ⊗ " + Pauli("XI").evolve(qc).to_label()[1])
    print("SWAP(I ⊗ X)SWAP† = " + Pauli("IX").evolve(qc).to_label()[0] + " ⊗ " + Pauli("IX").evolve(qc).to_label()[1])
    print("SWAP(Z ⊗ I)SWAP† = " + Pauli("ZI").evolve(qc).to_label()[0] + " ⊗ " + Pauli("ZI").evolve(qc).to_label()[1])
    print("SWAP(I ⊗ Z)SWAP† = " + Pauli("IZ").evolve(qc).to_label()[0] + " ⊗ " + Pauli("IZ").evolve(qc).to_label()[1])
    print("")
    
swap_gate()