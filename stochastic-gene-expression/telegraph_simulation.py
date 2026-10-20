"""Gillespie simulation of the telegraph model of gene expression.

Compares slow vs. fast promoter switching with the same mean mRNA level,
and checks the simulated Fano factor against the analytical prediction.
"""
import numpy as np
import matplotlib.pyplot as plt


def telegraph(k_on, k_off, k_tx, gamma, t_max, rng):
    """Simulate one cell. Returns (gene state, mRNA count) at time t_max."""
    t, g, m = 0.0, 0, 0
    while True:
        rates = np.array([k_on * (1 - g),   # OFF -> ON
                          k_off * g,         # ON -> OFF
                          k_tx * g,          # transcription (only when ON)
                          gamma * m])        # mRNA degradation
        total = rates.sum()
        dt = rng.exponential(1 / total)
        if t + dt > t_max:
            return g, m
        t += dt
        j = rng.choice(4, p=rates / total)  # which reaction fires
        if j == 0:
            g = 1
        elif j == 1:
            g = 0
        elif j == 2:
            m += 1
        else:
            m -= 1


def theory(k_on, k_off, k_tx, gamma):
    """Analytical steady-state mean and Fano factor of the telegraph model."""
    mean = (k_tx / gamma) * k_on / (k_on + k_off)
    fano = 1 + k_tx * k_off / ((k_on + k_off) * (k_on + k_off + gamma))
    return mean, fano


def main(n_cells=1000, seed=0):
    rng = np.random.default_rng(seed)
    params = {"slow switching": (0.05, 0.05, 50, 1),
              "fast switching": (5, 5, 50, 1)}
    results = {}
    for name, p in params.items():
        x = np.array([telegraph(*p, t_max=100, rng=rng)[1] for _ in range(n_cells)])
        results[name] = x
        mean_th, fano_th = theory(*p)
        print(f"{name}: mean={x.mean():.1f} (theory {mean_th:.1f}), "
              f"Fano={x.var() / x.mean():.1f} (theory {fano_th:.2f}), "
              f"fraction near 0={np.mean(x < 5):.2f}")

    for name, x in results.items():
        plt.hist(x, bins=40, alpha=0.6, label=name)
    plt.xlabel("mRNA count per cell")
    plt.ylabel("Number of cells")
    plt.title("Telegraph model: slow vs. fast promoter switching")
    plt.legend()
    plt.savefig("histogram.png", dpi=150, bbox_inches="tight")
    print("Saved histogram.png")


if __name__ == "__main__":
    main()
