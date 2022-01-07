import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def main():
    # arr = np.loadtxt("/home/jlb/times.txt")
    # tot = arr.sum(axis=1)
    # # print("\n".join(f"{t:.3f}" for t in tot))
    # res = np.reshape(tot[:90], (5, 18))
    # psi = res / res[0]
    # # print("\n".join(f"{t:.3f}" for t in np.reshape(psi, 90)))
    # last = tot[90:] / tot[:12]
    # # print("\n".join(f"{t:.3f}" for t in last))
    # p = np.concatenate([[1] * 18, [4] * 18, [9] * 18, [16] * 18, [28] * 18, [28] * 12])
    # psi = np.loadtxt("/home/jlb/psi.txt")
    # F = (1 / psi - 1 / p) / (1 - 1 / p)
    # # print("\n".join(f"{t:.3f}" for t in F))
    params = np.loadtxt("/home/jlb/p S P.txt")
    ys = np.loadtxt("/home/jlb/psi F.txt")
    data = np.concatenate((params, ys), axis=1)
    df = pd.DataFrame(data, columns=["ps", "S", "P", "psi", "F"])
    df = df.sort_values(by=["S", "P", "ps"])
    reshaped = np.reshape(df.to_numpy(), (18, 4, 5))

    num_colours = 18
    colors = plt.cm.Spectral(np.linspace(0, 1, num_colours))

    fig: plt.Figure
    ax1: plt.Axes
    ax2: plt.Axes
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8))
    ax1.set_prop_cycle('color', colors)
    ax2.set_prop_cycle('color', colors)

    ax1.set_xlabel("Number of threads")
    ax1.set_ylabel("Speedup $\\psi$")
    ax2.set_xlabel("Number of threads")
    ax2.set_ylabel("Karp-Flatt metric $F$")

    for group in reshaped:
        S = group[0, 1]
        P = group[0, 2]
        label = f"$S = {S}, P = {P}$"
        ps = group[:, 0]
        psis = group[:, 3]
        fs = group[:, 4]
        ax1.plot(ps, psis, label=label)
        ax2.plot(ps, fs, label=label)

    ax1.legend()
    ax2.legend()
    fig.tight_layout()
    fig.show()


if __name__ == "__main__":
    main()
