import numpy as np
import matplotlib.pyplot as plt


def main():
    sizes = [20, 25, 30, 35, 40]
    times = [0.0161, 0.0753, 0.2271, 0.5283, 1.1835]
    sixth = [t**(1/6) for t in times]
    fit = np.poly1d(np.polyfit(sizes, sixth, 1))
    fig: plt.Figure
    ax: plt.Axes
    fig, ax = plt.subplots(1, figsize=(6, 4))
    ax.set_xlabel("Space Station Size $S$")
    ax.set_ylabel(r"$(t_{\rm search})^{\frac{1}{6}}$ / ${\rm s}^{\frac{1}{6}}$")
    ax.plot(sizes, sixth, "k+:", label="Experimental data", )
    x1, x2 = ax.get_xlim()
    xs = np.linspace(x1, x2)
    ax.plot(
        xs,
        fit(xs),
        ls="--",
        c="0.6",
        label=f"Line of best fit\n$(t_{{\\rm search}})^{{\\frac{{1}}{{6}}}} = {fit.coeffs[0]:.3f}S {fit.coeffs[1]:+.3f}$"
    )
    ax.legend()
    fig.tight_layout()
    fig.show()


if __name__ == "__main__":
    main()
