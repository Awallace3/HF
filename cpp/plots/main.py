import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def hf_timings1(datapath="../outputs/t1.out", fn="images/t1.png", title="Water") -> None:
    """
    hf_timings1 plots t1 and t3 timings for speed up and efficiency
    """
    # Read datapath
    with open(datapath, "r") as f:
        d = f.readlines()
    threads = []
    times = []
    for l in d:
        if "Time (USR)" in l:
            times.append(float(l.split()[-1]))
        if "Num Threads: " in l:
            threads.append(int(l.split()[-1]))
    threads.pop(0)
    serial_time = times.pop(0)
    speedup = [serial_time / t for t in times]
    efficiency = [serial_time / (t * n) for t, n in zip(times, threads)]
    print(threads)
    print(times)
    print(speedup)
    print(efficiency)
    c1 = "tab:blue"
    c2 = "tab:orange"
    fix, ax1 = plt.subplots(dpi=300)
    ax1.plot(threads, speedup, "-o", label="Speed Up", color=c1)
    ax1.tick_params(axis="y", labelcolor=c1)

    ax2 = ax1.twinx()
    ax2.plot(threads, efficiency, "-o", label="Efficiency", color=c2)
    ax2.tick_params(axis="y", labelcolor=c2)
    ax1.set_xlabel("Number of Threads")
    ax1.set_ylabel("Speed Up", color=c1)
    ax2.set_ylabel("Efficiency", color=c2)
    plt.title(f"Speed Up and Efficiency: {title}")
    plt.savefig(fn)
    return


def fig2():
    # Read data
    df = pd.read_csv("data/fig2.csv", header=None, names=["x", "y"])
    df2 = pd.read_csv("data/fig2_hive.csv")
    # print("Hive1 Results")
    # df2['ys'] = df2.apply(lambda row: row['matmul'] / row['schulz'], axis=1)
    df2['ys'] = df2.apply(lambda row: row['matmul'] / row['ortho'], axis=1)

    # Plot
    fig, ax = plt.subplots(dpi=300)
    c1 = "tab:blue"
    c2 = "tab:orange"
    ax.plot(df["x"], df["y"], "o", color=c1)
    ax.set_xlabel("Matrix Size (n)")
    ax.tick_params(axis="y", labelcolor=c1)
    ax.set_ylabel("MatMul : Diagonalization", color=c1)

    ax2 = ax.twinx()
    ax2.plot(df2['n'], df2['ys'], "o", label="Hive Results", color=c2)
    ax2.set_ylabel("MatMul : Schulz Iteration", color=c2)
    ax2.set_ylim(0, 0.040)
    ax2.tick_params(axis="y", labelcolor=c2)

    fig.tight_layout()
    plt.savefig("plots/fig2_ortho.png")
    return


def fig7() -> None:
    df = pd.read_csv("data/t1/all.csv")
    # df['Total'] = df.apply(lambda row: row['ortho'] + row['allreduce'] + row[
    #     'host_device'] + row['mem_init'] + row['schulz'],
    #                        axis=1)

    fig, ax = plt.subplots(dpi=300)
    # ax.plot(df["n"], df["Total"], "-o", label="Total")
    ax.plot(df["n"], df["ortho"], "-o", label="Total")
    ax.plot(df["n"], df["allreduce"], "-o", label="AllReduce")
    ax.plot(df["n"], df["host_device"], "-o", label="Host to Device")
    ax.plot(df["n"], df["mem_init"], "-o", label="Memory Initialization")
    ax.plot(df["n"], df["schulz"], "-o", label="Schulz Iterations")
    ax.plot(df["n"], df["MatMul"], "-o", label="MatMul")
    ax.set_xlabel("MPI Tasks")
    ax.set_ylabel("Time (s)")
    ax.set_yscale("log")
    plt.legend()
    plt.savefig("plots/fig7.png")

    # Bar chart
    fig, ax = plt.subplots(dpi=300)
    df.drop(columns=["n"], inplace=True)
    df.drop(columns=["gpus"], inplace=True)
    labels = df.columns.values
    print(labels)
    labels = [
        'Total',
        'AllReduce',
        'H to D',
        'Mem Init',
        'Schulz',
        'A^TA',
        "MatMul",
    ]
    colors = [
        'tab:blue',
        'tab:orange',
        'tab:green',
        'tab:red',
        'tab:purple',
        'tab:brown',
        'tab:pink',
    ]
    print(labels)
    values = df.iloc[0].values
    # ax.bar(labels, values, colors=colors)
    ax.bar(labels, values, color=colors)
    plt.xlabel("Tasks")
    plt.ylabel("Time (s)")
    ax.set_yscale("log")
    plt.savefig("plots/fig7_bar.png")
    return


def main():
    # fig2()
    # fig7()
    hf_timings1(datapath="../outputs/t1.out", fn="images/t1.png", title="Water HF/STO-3G")
    hf_timings1(datapath="../outputs/t3.out", fn="images/t3.png", title="Ethene HF/aug-cc-pVDZ")
    return


if __name__ == "__main__":
    main()
