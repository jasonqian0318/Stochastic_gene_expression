# Stochastic Gene Expression Simulation

Gillespie (stochastic simulation algorithm) implementation of the **telegraph model** of gene expression, where a gene switches between ON and OFF states and produces mRNA only when ON.

## What it shows

Two parameter sets with the **same mean** mRNA level (about 25 per cell), 1,000 simulated cells each:

- **Slow switching** (k_on = k_off = 0.05): bursty, bimodal distribution, high Fano factor.
- **Fast switching** (k_on = k_off = 5): unimodal distribution close to a simple birth–death process, low Fano factor.

The simulated mean and Fano factor (variance / mean) are compared with the analytical predictions:

- mean = (k_tx / γ) · k_on / (k_on + k_off)
- Fano = 1 + k_tx · k_off / ((k_on + k_off)(k_on + k_off + γ))

## How to run

```bash
pip install -r requirements.txt
python telegraph_simulation.py
```

The script prints the simulated vs. theoretical statistics and saves `histogram.png`.

## Model

| Reaction | Rate |
|---|---|
| OFF → ON | k_on |
| ON → OFF | k_off |
| ON → ON + mRNA | k_tx |
| mRNA → ∅ | γ · m |
