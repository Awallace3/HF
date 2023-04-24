import psi4
import numpy as np
import pandas as pd
from qm_tools_aw import tools
import qcelemental as qcel


def psi4_compute(mol, outdata="t2"):
    psi4.set_memory('2 GB')
    psi4.set_num_threads(10)
    psi4.core.set_output_file(f'{outdata}/output.dat', False)
    psi4.set_options({"basis": "aug-cc-pvdz"})
    wfn = psi4.core.Wavefunction.build(mol,
                                       psi4.core.get_global_option("basis"))
    mints = psi4.core.MintsHelper(wfn.basisset())
    S = np.asarray(mints.ao_overlap())
    np.savetxt(f"{outdata}/S.csv", S, delimiter=" ")
    # print(f"{S = }")
    T = np.asarray(mints.ao_potential())
    np.savetxt(f"{outdata}/T.csv", T, delimiter=" ")
    # print(f"{T = }")
    V = np.asarray(mints.ao_kinetic())
    np.savetxt(f"{outdata}/V.csv", V, delimiter=" ")
    # print(f"{V = }")
    I = np.asarray(mints.ao_eri())
    nbf = len(I)
    print(f"{nbf = }")

    with open(f"{outdata}/eri.csv", 'w') as f:
        for i in range(nbf):
            for j in range(i + 1):
                for k in range(i + 1):
                    for l in range(k + 1):
                        # print(i, j, k, l, I[i,j,k,l])
                        line = f"{i} {j} {k} {l} {I[i,j,k,l]}\n"
                        # line = f"{i} {j} {k} {l}\n"
                        f.write(line)

    e = psi4.energy("HF/aug-cc-pvdz")
    print(f"{e =}")
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
            return tools.print_cartesians_pos_carts(mmA[:,0], mmA[:,1:])

def main():
    # mol = qcel.models.Molecule.from_file("t1/t1.xyz")
    d = find_geoms(6)
    # fn = "t1/t1.xyz"
    # with open(fn, "r") as f:
    #     d = "".join(f.readlines()[2:])
    #     print(d)
    print(d)
    mol = psi4.geometry(d)
    psi4_compute(mol)
    return


if __name__ == "__main__":
    main()
