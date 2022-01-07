import numpy as np
import matplotlib.pyplot as plt


def main():
    sizes = [20, 25, 30, 35, 40]
    times_1 = np.array([0.0161, 0.0753, 0.2271, 0.5283, 1.1835]) ** (1/6)
    times_5 = np.array([0.0005, 0.0017, 0.0052, 0.0075, 0.0225]) ** (1/6)
    times_10 = np.array([0.0002, 0.0004, 0.0010, 0.0020, 0.0032]) ** (1/6)
    times_15 = np.array([0.0001, 0.0003, 0.0005, 0.0012, 0.0016]) ** (1/6)
    times_20 = np.array([0.0001, 0.0002, 0.0003, 0.0007, 0.0014]) ** (1/6)
    fig: plt.Figure
    ax1: plt.Axes
    ax2: plt.Axes
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    ax1.set_xlabel("Space Station Size $S$")
    ax1.set_ylabel(r"$(t_{\rm search})^{\frac{1}{6}}$ / ${\rm s}^{\frac{1}{6}}$")
    ax1.plot(sizes, times_5, "+-", c="0.6", label="${\\tt n\\_cells} = 5$")
    ax1.plot(sizes, times_10, "+:", c="0.4", label="${\\tt n\\_cells} = 10$")
    ax1.plot(sizes, times_15, "+--", c="0.2", label="${\\tt n\\_cells} = 15$")
    ax1.plot(sizes, times_20, "+-.", c="0.0", label="${\\tt n\\_cells} = 20$")
    ax1.legend()
    ax2.set_xlabel("Space Station Size $S$")
    ax2.set_ylabel(r"$(t_{\rm search})^{\frac{1}{6}}$ / ${\rm s}^{\frac{1}{6}}$")
    ax2.plot(sizes, times_1, "+-", c="0.8", label="${\\tt n\\_cells} = 1$")
    ax2.plot(sizes, times_5, "+-", c="0.6", label="${\\tt n\\_cells} = 5$")
    ax2.plot(sizes, times_10, "+:", c="0.4", label="${\\tt n\\_cells} = 10$")
    ax2.plot(sizes, times_15, "+--", c="0.2", label="${\\tt n\\_cells} = 15$")
    ax2.plot(sizes, times_20, "+-.", c="0.0", label="${\\tt n\\_cells} = 20$")
    ax2.legend()
    fig.tight_layout()
    fig.show()


if __name__ == "__main__":
    main()
