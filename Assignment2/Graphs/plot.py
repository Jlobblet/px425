import pandas as pd
import matplotlib.pyplot as plt


def main():
    plt.style.use("grayscale")
    plt.rcParams["savefig.facecolor"] = "none"

    df = pd.read_csv("data.csv")
    compilers = ["gcc", "clang", "icc"]
    styles = ["-", "--", "-."]
    optimisation_levels = ["0", "1", "2", "3"]

    for o in optimisation_levels:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.grid(ls=":")
        for compiler, style in zip(compilers, styles):
            ax.plot(df[f"{compiler}{o}"], label=compiler, ls=style)

        ax.set_xlabel("Index of Change")
        ax.set_ylabel("Execution Time / $s$")
        ax.legend(loc="upper right")
        fig.suptitle(f"Comparison of execution time on -O{o}")
        fig.tight_layout()
        fig.savefig(f"{o}.png")


if __name__ == "__main__":
    main()
