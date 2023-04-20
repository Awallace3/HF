import numpy as np
import math
import sys
import re
import pandas as pd


def clean_line(line):
    for _ in range(5):
        line = line.replace("  ", " ")
    line = line.split(" ")
    line = line[-1].rstrip()
    return float(line)


def clean_line_last(line):
    for _ in range(5):
        line = line.replace("  ", " ")
    line = line.split(" ")
    length = line[1].rstrip()
    line = line[-1].rstrip()

    return float(line), int(length)


def eigen(A):
    # since np.linalg.eig doesn't ensure sorting occurs, the following code has been added from a stackoverflow post
    # https://stackoverflow.com/questions/8092920/sort-eigenvalues-and-associated-eigenvectors-after-using-numpy-linalg-eig-in-pyt
    eigenValues, eigenVectors = np.linalg.eig(A)
    idx = np.argsort(eigenValues)
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:, idx]
    return (eigenValues, eigenVectors)


def parse_input(path):
    with open(path, "r") as fp:
        v = fp.readlines()
        val, size = clean_line_last(v[-1])
        vs = np.zeros((size, size))
    for i in range(len(v)):
        l = v[i].rstrip()
        while "  " in l:
            l = l.replace("  ", " ")
        if l[0] == " ":
            l = l[1:].split(" ")  # first spot is whitespace
        else:
            l = l.split(" ")
        l[0] = int(l[0]) - 1
        l[1] = int(l[1]) - 1
        l[2] = float(l[2])
        vs[l[0], l[1]] = l[2]
        vs[l[1], l[0]] = l[2]

    return vs, size


def read_inputs(input_dir):
    with open("%s/enuc.dat" % input_dir, "r") as fp:
        enuc = float(fp.read().rstrip())
    s, l1 = parse_input("%s/s.dat" % input_dir)
    t, l2 = parse_input("%s/t.dat" % input_dir)
    v, l3 = parse_input("%s/v.dat" % input_dir)
    if l2 != l3:
        print("T and V have different lengths. Fix your input files.")
        sys.exit()

    return enuc, s, t, v, l2


def print_2d_array(array, size):
    df = pd.DataFrame(array, columns=["%s" % i for i in range(1, size + 1)])
    df.index += 1
    pd.set_option("display.float_format", lambda x: "%.7f" % x)
    print(df)


def orthogonalization_of_basis_set(s):
    # snippet from pentavalentcarbon, https://chemistry.stackexchange.com/questions/85484/how-to-do-lowdin-symmetric-orthonormalisation/85485
    lam_s, l_s = np.linalg.eig(s)
    lam_s = lam_s * np.eye(len(lam_s))
    lam_sqrt_inv = np.sqrt(np.linalg.inv(lam_s))
    symm_orthog = np.dot(l_s, np.dot(lam_sqrt_inv, l_s.T))
    # end of snippet
    return symm_orthog


def initial_Fock_matrix(s_orth, H_core):
    symm_orth_t = np.transpose(s_orth)
    s_h = np.matmul(symm_orth_t, H_core)
    F_0 = np.matmul(s_h, s_orth)
    return F_0


def C_prime_gen(F):
    epsilon, C_prime = np.linalg.eig(F)
    idx = np.argsort(epsilon)
    epsilon = epsilon[idx]
    C_prime = C_prime[:, idx]
    return C_prime


def Fock_Matrix(s_orth, F):
    symm_orth_t = np.transpose(s_orth)
    s_h = np.matmul(symm_orth_t, F)
    F_prime = np.matmul(s_h, s_orth)
    return F_prime


def initial_Density_matrix_einsum(s_orth, H_core):
    # https://ajcr.net/Basic-guide-to-einsum/
    s_orth_t = np.transpose(s_orth)
    s_h = np.einsum("ia,ab->ib", s_orth_t, H_core)
    F_0 = np.einsum("ib,bj->ij", s_h, s_orth)
    return F_0


def Density_Matrix(C_0, e):
    occ = math.ceil(e / 2)
    # D_einsum = np.einsum("ij, ik->jk", C_0, C_0)
    D = np.zeros(np.shape(C_0))
    # D_mu,nu =
    for m in range(occ):
        for mu in range(len(D)):
            for nu in range(len(D)):
                D[mu, nu] += C_0[mu, m] * C_0[nu, m]
    """
    #
    for mu in range(len(D)):
        for nu in range(len(D)):
            for m in range(occ):
                tmp += (C_0[mu, m] * C_0[nu, m])
            D[mu, nu] = tmp
    """
    return D


def initial_E_elec(D_0, H):
    E_0_elec = 0
    for mu in range(len(H)):
        for nu in range(len(H)):
            E_0_elec += D_0[mu, nu] * (2 * H[mu, nu])
    return E_0_elec


def E_elec_gen(D_0, H, F):
    E_0_elec = 0
    for mu in range(len(H)):
        for nu in range(len(H)):
            E_0_elec += D_0[mu, nu] * (H[mu, nu] + F[mu, nu])
    return E_0_elec


def mu_nu_la_si_converter(mu, nu, la, si):
    mu_nu = mu * (mu + 1) / 2 + nu
    la_si = la * (la + 1) / 2 + si
    mu_nu_la_si = int(mu_nu * (mu_nu + 1) / 2 + la_si)
    return mu_nu_la_si


def Index(x, y):
    if x > y:
        return int((x * (x + 1) / 2) + y)
    else:
        return int((y * (y + 1) / 2) + x)


def read_two_e_repulsion(input_dir, size):
    # should be 4-D array but only giving into 1-d array

    eri = np.genfromtxt("%s/eri.dat" % input_dir)
    # M = ( 4 * size)
    M = (size) ** 2
    one_d_size = int(M * (M + 1) / 2)
    eri_vals = np.zeros((one_d_size))
    for row in eri:
        mu = int(row[0]) - 1
        nu = int(row[1]) - 1
        la = int(row[2]) - 1
        si = int(row[3]) - 1
        val = row[4]
        ij = Index(mu, nu)
        kl = Index(la, si)
        ijkl = Index(ij, kl)

        eri_vals[ijkl] = val

    # print('num zeros')
    # print(len(np.where(eri_vals==0)[0]))
    # print(one_d_size - 228)
    # print(eri_vals)
    return eri_vals


def new_Fock_Matrix(H, D, eri, size):
    e = 10
    occ = range(5)

    s = range(size)
    F = np.zeros(np.shape(H))
    for i in s:
        for j in s:
            for k in s:
                for l in s:
                    ij = Index(i, j)
                    kl = Index(k, l)
                    ijkl = Index(ij, kl)
                    ik = Index(i, k)
                    jl = Index(j, l)
                    ikjl = Index(ik, jl)

                    F[i, j] += D[k, l] * (2 * eri[ijkl] - eri[ikjl])

            # F[i, j] += H[i, j]
    # print(F)
    F = np.add(F, H)
    return F


def rms_Density_Matrix(D_new, D_old):
    s = range(len(D_old))
    val = 0
    for mu in s:
        for nu in s:
            val += math.pow((D_new[mu, nu] - D_old[mu, nu]), 2)
    return val ** (1 / 2)


def self_consistent_field_iteration(
    H, D_0, symm_orth, eri, e, enuc, E_0_elec, size, thres1, thres2
):

    iterations = 1
    D_old = D_0
    E_elec_old = E_0_elec
    rms_d = thres1 + 1
    delta_E = thres2 + 2
    E_elec_new = 0
    E_0_tot = E_elec_old + enuc

    print("Iter", "\tE(elec)", "\tE(tot)", "\tDelta(E)", "\tRMS(D)")
    print(iterations, E_elec_old, E_0_tot)
    not_converged = True
    while not_converged:
        F = new_Fock_Matrix(H, D_old, eri, size)
        # print("\nFock Matrix")
        # print_2d_array(F, size)

        F_prime = Fock_Matrix(symm_orth, F)
        C_prime = C_prime_gen(F_prime)
        C = np.matmul(symm_orth, C_prime)
        D_new = Density_Matrix(C, e)
        # print("\nD:")
        # print_2d_array(D_new, size)

        E_elec_new = E_elec_gen(D_new, H, F)
        rms_d = rms_Density_Matrix(D_old, D_new)
        delta_E = abs(E_elec_new - E_elec_old)
        E_tot = E_elec_new + enuc
        # print("\nIter\t\tE(elec)\t\tE(tot)\t\tDelta(E)\t\tRMS(D)")
        print(iterations, E_elec_new, E_tot, delta_E, rms_d)

        if (rms_d < thres1) and (delta_E < thres2):
            not_converged = False
        else:
            D_old = D_new
            E_elec_old = E_elec_new
            iterations += 1

    print("HF Energy", E_tot)
    print("Iterations", iterations)


def main():
    input_dir = "input/h2o/STO-3G"
    # input_dir = "input/h2o/DZ"
    e = 10
    thres1 = 1e-11
    thres2 = 1e-11

    enuc, s, t, v, size = read_inputs(input_dir)
    # Core Hamiltonian matrix, H
    H = np.add(t, v)
    print("\nH")
    print_2d_array(H, size)

    symm_orth = orthogonalization_of_basis_set(s)
    print("\nS^(-1/2)")
    print_2d_array(symm_orth, size)

    F_0 = initial_Fock_matrix(symm_orth, H)
    print("\nF_0")
    print_2d_array(F_0, size)

    # III. 2
    C_0_prime = C_prime_gen(F_0)

    # print("C_0_prime")
    print("\nC_0_prime:")
    print_2d_array(C_0_prime, size)

    # III. 3
    C_0 = np.matmul(symm_orth, C_0_prime)

    print("\nC_0")
    print_2d_array(
        C_0, size
    )  # brent said sign doesn't matter here so try next
    # III. 4
    D_0 = Density_Matrix(C_0, e)
    print("\nD_0")
    print_2d_array(D_0, size)

    # III. 5
    E_0_elec = initial_E_elec(D_0, H)
    print("\nE_0_elec", E_0_elec)
    SCF_0 = E_0_elec + enuc
    print("SCF_0", SCF_0, "\n")

    print("Size", size)
    eri = read_two_e_repulsion(input_dir, size)

    self_consistent_field_iteration(
        H, D_0, symm_orth, eri, e, enuc, E_0_elec, size, thres1, thres2
    )


if __name__ == "__main__":
    main()
