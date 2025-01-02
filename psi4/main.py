import psi4
import numpy as np
import os


def psi4_compute(
    mol,
    outdata="t2",
    level_of_theory="HF/aug-cc-pvdz",
    mol_units="units=angstrom",
) -> None:
    if not os.path.exists(outdata):
        os.makedirs(outdata)
    with open(f"{outdata}/geom.xyz", "w") as f:
        n = mol.count("\n")
        f.write(f"{n - 1}\n")
        f.write(mol)

    mol = psi4.geometry(mol)
    psi4.set_memory("60 GB")
    psi4.set_num_threads(20)
    psi4.core.set_output_file(f"{outdata}/output.dat", False)
    method, basis = level_of_theory.split("/")
    psi4.set_options({"basis": basis})
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("basis"))
    mints = psi4.core.MintsHelper(wfn.basisset())
    S = np.asarray(mints.ao_overlap())
    np.savetxt(f"{outdata}/S.csv", S, delimiter=" ", fmt="%1.7f")
    T = np.asarray(mints.ao_potential())
    np.savetxt(f"{outdata}/T.csv", T, delimiter=" ", fmt="%1.7f")
    V = np.asarray(mints.ao_kinetic())
    np.savetxt(f"{outdata}/V.csv", V, delimiter=" ", fmt="%1.7f")
    eris = np.asarray(mints.ao_eri())
    nbf = len(eris)
    print(f"{nbf=}")
    enuc = mol.nuclear_repulsion_energy()
    print(f"{enuc=}")
    with open(f"{outdata}/enuc.csv", "w") as f:
        f.write(f"{enuc}")

    with open(f"{outdata}/eri.csv", "w") as f:
        for i in range(nbf):
            for k in range(i + 1):
                for j in range(i + 1):
                    for l in range(k + 1):
                        line = f"{i} {j} {k} {l} {eris[i, j, k, l]:.15f}\n"
                        f.write(line)

    e = psi4.energy(level_of_theory, molecule=mol)
    print(f"Target Energy ({level_of_theory}): {e}")
    with open(f"{outdata}/final_energy.dat", "w") as f:
        f.write(f"{e}\n")
    return


def benzene():
    return """
6       1.5000000000    -1.8000000000   -1.3915000000
6       2.7050743494    -1.8000000000   -0.6957500000
6       2.7050743494    -1.8000000000   0.6957500000
6       1.5000000000    -1.8000000000   1.3915000000
6       0.2949256506    -1.8000000000   0.6957500000
6       0.2949256506    -1.8000000000   -0.6957500000
1       1.5000000000    -1.8000000000   -2.4715000000
1       3.6403817855    -1.8000000000   -1.2357500000
1       3.6403817855    -1.8000000000   1.2357500000
1       1.5000000000    -1.8000000000   2.4715000000
1       -0.6403817855   -1.8000000000   1.2357500000
1       -0.6403817855   -1.8000000000   -1.2357500000
"""


def t1_example():
    return """
8    0.000000000000    0.000000000000    -0.071151380605
1    0.000000000000    0.757939245855    0.564612021746
1    0.000000000000    -0.757939245855    0.564612021746
    """


def main():
    d = t1_example()
    psi4_compute(d, outdata="../cpp/data/t1", level_of_theory="HF/STO-3G")
    return


if __name__ == "__main__":
    main()
