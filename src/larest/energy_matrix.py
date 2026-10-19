"""Energy descriptors of ring-opening model reactions, collected in energy matrices.

For a monomer M, polymer model compounds P_n of degree n (linear for ROR,
cyclic for RER) and, for ROR, an initiator XY, the following descriptors are
computed for each thermodynamic quantity X (E, H, or G):

* ``X^n_n = X(P_n) - n X(M) - X(XY)`` — energy of the model reaction forming
  P_n (no initiator term for RER).
* ``X^n_0 = X^n_n / n`` — the reaction energy normalised per monomer.
* ``IX^n_m = (X(P_n) - X(P_m)) / (n - m) - X(M)`` for ``0 < m < n`` — the
  incremental polymerisation energy per monomer for inserting ``n - m``
  monomers into a chain of degree m.  For RER, the cyclic "polymer" of degree 1
  is the monomer itself.

These are arranged in a lower-triangular matrix with rows n and columns m:
column 0 holds ``X^n_0``, the diagonal holds ``X^n_n`` and the off-diagonal
entries hold ``IX^n_m``.  A final row gives ``X^∞_0``, extrapolated to infinite
chain length by linear regression of ``X^n_0`` against 1/n, together with the
coefficient of determination R² of the fit.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

R2_WARNING_THRESHOLD = 0.8
INF_LABEL = "inf"
R2_COLUMN = "R2"


def extrapolate_to_infinite_length(
    lengths: Sequence[float],
    values: Sequence[float],
) -> tuple[float, float | None]:
    """Extrapolate per-monomer values to infinite chain length.

    Fits ``value = a / n + b`` by linear regression of *values* against
    ``1 / n`` and returns the intercept ``b``.

    Parameters
    ----------
    lengths : list[float]
        Chain lengths n.
    values : list[float]
        Values at each chain length (NaN values are ignored).

    Returns
    -------
    tuple[float, float | None]
        The extrapolated value (NaN if there are fewer than two valid points)
        and the R² of the fit (``None`` if fewer
        than three points were fitted, since two points always fit exactly).
    """
    x = np.array(lengths, dtype=np.float64)
    y = np.array(values, dtype=np.float64)
    valid = ~np.isnan(x) & ~np.isnan(y)
    x, y = x[valid], y[valid]
    if len(x) < 2:
        # no extrapolation is possible from a single chain length
        return np.nan, None

    slope, intercept = np.polyfit(1.0 / x, y, 1)
    if len(x) == 2:
        # a line through two points always fits exactly; R² is meaningless
        return float(intercept), None
    residuals = y - (slope / x + intercept)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(residuals**2)) / ss_tot if ss_tot > 0 else 1.0
    return float(intercept), r2


def compute_energy_matrix(
    monomer: float | None,
    polymers: Mapping[int, float | None],
    initiator: float | None,
    reaction_type: str,
) -> pd.DataFrame:
    """Compute the energy matrix for one thermodynamic quantity.

    Parameters
    ----------
    monomer : float | None
        Value for the monomer (J/mol).
    polymers : dict[int, float | None]
        Values for each polymer model compound, keyed by degree n (J/mol).
    initiator : float | None
        Value for the initiator (J/mol); ignored for RER.
    reaction_type : str
        ``"ROR"`` or ``"RER"``.

    Returns
    -------
    pd.DataFrame
        Matrix in kJ/mol with one row per degree n (plus an ``"inf"`` row) and
        columns ``0 .. max(n)`` for m, plus an ``"R2"`` column that is set on the
        ``"inf"`` row only.  Entries that cannot be computed are NaN.
    """
    lengths = sorted(polymers)
    max_length = max(lengths, default=0)
    columns = [*range(max_length + 1), R2_COLUMN]
    matrix = pd.DataFrame(
        np.nan,
        index=pd.Index([*lengths, INF_LABEL], name="n", dtype=object),
        columns=pd.Index(columns, name="m", dtype=object),
    )

    def value(x: float | None) -> float:
        return np.nan if x is None else x

    m_value = value(monomer)
    i_value = value(initiator) if reaction_type == "ROR" else 0.0
    chain_values: dict[int, float] = {n: value(x) for n, x in polymers.items()}
    if reaction_type == "RER":
        chain_values.setdefault(1, m_value)

    for n in lengths:
        reaction = chain_values[n] - n * m_value - i_value
        matrix.loc[n, n] = reaction / 1000
        matrix.loc[n, 0] = reaction / n / 1000
        for m in range(1, n):
            if m in chain_values:
                increment = (chain_values[n] - chain_values[m]) / (n - m) - m_value
                matrix.loc[n, m] = increment / 1000

    extrapolated, r2 = extrapolate_to_infinite_length(
        lengths,
        [matrix.loc[n, 0] for n in lengths],
    )
    matrix.loc[INF_LABEL, 0] = extrapolated
    matrix.loc[INF_LABEL, R2_COLUMN] = np.nan if r2 is None else r2
    return matrix


def compute_energy_matrices(
    monomer: dict[str, float | None],
    polymers: dict[int, dict[str, float | None]],
    initiator: dict[str, float | None] | None,
    reaction_type: str,
    params: list[str],
) -> dict[str, pd.DataFrame]:
    """Compute energy matrices for several thermodynamic quantities.

    Parameters
    ----------
    monomer : dict[str, float | None]
        Monomer parameters for one pipeline section.
    polymers : dict[int, dict[str, float | None]]
        Polymer parameters for the same section, keyed by degree n.
    initiator : dict[str, float | None] | None
        Initiator parameters for the same section (ROR only).
    reaction_type : str
        ``"ROR"`` or ``"RER"``.
    params : list[str]
        Thermodynamic quantities to compute matrices for, e.g. ``["E", "H", "G"]``.

    Returns
    -------
    dict[str, pd.DataFrame]
        Mapping of quantity to its energy matrix (see
        :func:`compute_energy_matrix`).
    """
    return {
        param: compute_energy_matrix(
            monomer=monomer.get(param),
            polymers={n: values.get(param) for n, values in polymers.items()},
            initiator=initiator.get(param) if initiator else None,
            reaction_type=reaction_type,
        )
        for param in params
    }


def format_energy_matrix(matrix: pd.DataFrame, indent: str = "") -> str:
    """Format an energy matrix as a human-readable table string.

    Parameters
    ----------
    matrix : pd.DataFrame
        Matrix as returned by :func:`compute_energy_matrix`.
    indent : str
        String prepended to every line.

    Returns
    -------
    str
        Table with rows n, columns m, and R² shown after the extrapolated value.
    """
    display = matrix.drop(columns=R2_COLUMN).copy()
    display.index = [("∞" if n == INF_LABEL else str(n)) for n in display.index]
    display.columns = [f"m={m}" for m in display.columns]
    lines = display.to_string(float_format="{:.2f}".format, na_rep="").split("\n")
    r2 = matrix.loc[INF_LABEL, R2_COLUMN]
    if not np.isnan(r2):
        lines[-1] += f"   (R² = {r2:.3f})"
    lines.insert(len(lines) - 1, "-" * len(lines[0]))
    return "\n".join(f"{indent}{line}" for line in lines)


def low_r2_warning(matrix: pd.DataFrame) -> bool:
    """Return ``True`` if the extrapolation R² is below the warning threshold."""
    r2 = matrix.loc[INF_LABEL, R2_COLUMN]
    return bool(not np.isnan(r2) and r2 < R2_WARNING_THRESHOLD)
