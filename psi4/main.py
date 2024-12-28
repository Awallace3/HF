import psi4
import numpy as np
import pandas as pd
from qm_tools_aw import tools
import qcelemental as qcel
import os


def psi4_compute(mol, outdata="t2", level_of_theory='HF/aug-cc-pvdz', mol_units="units=angstrom") -> None:
    if not os.path.exists(outdata):
        os.makedirs(outdata)
    with open(f"{outdata}/geom.xyz", 'w') as f:
        n = mol.count('\n')
        f.write(f"{n - 1}\n")
        f.write(mol)

    mol = psi4.geometry(mol)
    psi4.set_memory('60 GB')
    psi4.set_num_threads(20)
    psi4.core.set_output_file(f'{outdata}/output.dat', False)
    method, basis = level_of_theory.split("/")
    psi4.set_options({"basis": basis})
    wfn = psi4.core.Wavefunction.build(mol,
                                       psi4.core.get_global_option("basis"))
    mints = psi4.core.MintsHelper(wfn.basisset())
    S = np.asarray(mints.ao_overlap())
    np.savetxt(f"{outdata}/S.csv", S, delimiter=" ", fmt='%1.7f')
    T = np.asarray(mints.ao_potential())
    np.savetxt(f"{outdata}/T.csv", T, delimiter=" ", fmt='%1.7f')
    V = np.asarray(mints.ao_kinetic())
    np.savetxt(f"{outdata}/V.csv", V, delimiter=" ", fmt='%1.7f')

    I = np.asarray(mints.ao_eri())
    nbf = len(I)
    print(f"{nbf = }")

    enuc = mol.nuclear_repulsion_energy()
    print(f"{enuc = }")
    with open(f"{outdata}/enuc.csv", 'w') as f:
        f.write(f"{enuc}")


    with open(f"{outdata}/eri.csv", 'w') as f:
        for i in range(nbf):
            for k in range(i + 1):
                for j in range(i + 1):
                    for l in range(k + 1):
                        line = f"{i} {j} {k} {l} {I[i,j,k,l]:.15f}\n"
                        f.write(line)

    e = psi4.energy(level_of_theory, molecule=mol)
    print(f"Target Energy ({level_of_theory}): {e}")
    with open(f"{outdata}/final_energy.dat", 'w') as f:
        f.write(f"{e}\n")
    return

def find_geoms(size=10) -> None:
    """
    find_geoms
    """
    mol = psi4.geometry("""
0 1
8    0.000000000000    0.000000000000    -0.071151380605
1    0.000000000000    0.757939245855    0.564612021746
1    0.000000000000    -0.757939245855    0.564612021746
    """)
    df = pd.read_pickle("schr.pkl")
    print(df)
    for n, r in df.iterrows():
        if len(r['monAs']) == size:
            mmA = r['Geometry'][r['monAs'], :]
            # tools.print_cartesians_pos_carts(mmA[:,0], mmA[:,1:])
            return tools.print_cartesians_pos_carts(mmA[:,0], mmA[:,1:])

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

def main_example():
    return """
0 1
8   -0.702196054   -0.056060256   0.009942262
1   -1.022193224   0.846775782   -0.011488714
1   0.257521062   0.042121496   0.005218999
--
0 1
8   2.268880784   0.026340101   0.000508029
1   2.645502399   -0.412039965   0.766632411
1   2.641145101   -0.449872874   -0.744894473
units angstrom
"""

def t1_example():
    return """
8    0.000000000000    0.000000000000    -0.071151380605
1    0.000000000000    0.757939245855    0.564612021746
1    0.000000000000    -0.757939245855    0.564612021746
    """

def main():
    # d = find_geoms(6) # Ethene
    # return
    # d = benzene()
    # d = t1_example()
    # psi4_compute(d, outdata="../cpp/data/example1", level_of_theory='HF/cc-pvdz')
    d = t1_example()
    psi4_compute(d, outdata="../cpp/data/t1", level_of_theory='HF/STO-3G')
    return


if __name__ == "__main__":
    main()
