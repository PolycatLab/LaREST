# LaREST

**LaREST** (Lactone Ring-opening Energetics Sorting Tool) is a computational chemistry pipeline for predicting the thermodynamics of lactone ring-opening polymerization reactions.

Given a set of monomer SMILES strings and polymer chain lengths, it computes electronic energies (E), enthalpies (H), entropies (S), and free energies (G) — for the lowest free energy conformer and as Boltzmann averages — by chaining together four levels of theory:

| Stage | Tool | Description |
|---|---|---|
| 1 | RDKit + xTB | MMFF conformer generation, ranked by xTB free energy |
| 2 | CREST | Conformer/rotamer ensemble exploration (iMTD-GC) |
| 3 | CENSO + ORCA | DFT refinement of the ensemble (four sub-stages) |
| 4 | CREST entropy | Conformational entropy correction |

Two reaction types are supported:

- **RER** (Ring Expansion Reaction) — no initiator; cyclic oligomer
- **ROR** (Ring-Opening Reaction) — includes an initiator (water, alcohol, or an isodesmic ester/thioester/amide chosen automatically)

Final reaction thermodynamics are computed as:

$$\Delta X = \frac{X_\text{polymer} - n \cdot X_\text{monomer} - X_\text{initiator}}{n}$$

where *n* is the polymer chain length. Incremental polymerisation energies and
extrapolations to infinite chain length are collected in energy matrices (see
[Output](output.md#energy-matrices)).

```{toctree}
:maxdepth: 2
:caption: Getting started
:hidden:

installation
usage
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

configuration
pipeline
output
testing
```
