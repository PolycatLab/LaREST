"""Physical constants and shared lookup tables used throughout the LaREST pipeline.

Unit conversion factors convert values produced by external QM programs (xTB,
CENSO) into SI-consistent J/mol units before they are stored or compared.

SMARTS patterns in ``MONOMER_GROUPS`` and ``INITIATOR_GROUPS`` define the
ring-opening functional groups that LaREST recognises when building polymer
chains.
"""

HARTTREE_TO_JMOL = 2625499.63948
KCALMOL_TO_JMOL = 4184.0
CALMOL_TO_JMOL = 4.184
GAS_CONSTANT = 8.314462618  # J/mol/K
# Ring-opening functional groups, in priority order (highest first). When a
# monomer contains several breakable bonds, the group listed first wins. Each
# pattern maps the two atoms of the bond cleaved on ring opening: ``:1`` is the
# acyl (C=X) carbon, attacked by the initiator nucleophile, and ``:2`` is the
# ring atom it is bonded to, which becomes the new chain end. ``@`` restricts
# both bonds to ring bonds.
MONOMER_GROUPS: dict[str, str] = {
    "L": "[C,c;R]@[C;R:1](=[O;!R])@[O;R:2]",  # lactone
    "tnL": "[C,c;R]@[C;R:1](=[S;!R])@[O;R:2]",  # thionolactone
    "tL": "[C,c;R]@[C;R:1](=[O;!R])@[S;R:2]",  # thiolactone
    "dtL": "[C,c;R]@[C;R:1](=[S;!R])@[S;R:2]",  # dithiolactone
    "Lm": "[C,c;R]@[C;R:1](=[O;!R])@[N;R:2]",  # lactam
    "tLm": "[C,c;R]@[C;R:1](=[S;!R])@[N;R:2]",  # thiolactam
    # carbonates and related derivatives (not yet benchmarked)
    "CC": "[O;R]@[C;R:1](=[O;!R])@[O;R:2]",  # carbonate
    "CtC": "[S;R]@[C;R:1](=[O;!R])@[O;R:2]",  # thiocarbonate
    "CdtC": "[S;R]@[C;R:1](=[O;!R])@[S;R:2]",  # dithiocarbonate
    "CtnC": "[O;R]@[C;R:1](=[S;!R])@[O;R:2]",  # thionocarbonate
    "CX": "[S;R]@[C;R:1](=[S;!R])@[O;R:2]",  # xanthate
    "CtX": "[S;R]@[C;R:1](=[S;!R])@[S;R:2]",  # trithiocarbonate
    "oA": "[N;R]@[C;R:1](=[O;!R])@[O;R:2]",  # cyclic carbamate (urethane)
}
# Initiator groups, in priority order. ``:1`` is the nucleophilic atom that
# bonds to the monomer acyl carbon and ``:2`` is the atom (H or acyl C) that caps
# the ring heteroatom end of the chain. Explicit hydrogens are added before
# matching.
INITIATOR_GROUPS: dict[str, str] = {
    # water / alcohols (hydrolysis / alcoholysis); excludes carboxylic acids
    "R-OH": "[O;X2;H1,H2;!$(O[C,c]=[O,S]):1]-[#1:2]",
    # esters, thioesters and amides (isodesmic reactions), e.g. methyl acetate:
    # the alkoxy O caps the acyl end and the acyl C caps the ring heteroatom
    # end, e.g. CH3C(=O)-O-[chain]-C(=O)-O-CH3
    "R-COO-R": "[C;!R:2](=[O,S;!R])-[O;X2;!R:1]-[#6]",
    "R-COS-R": "[C;!R:2](=[O,S;!R])-[S;X2;!R:1]-[#6]",
    "R-CON-R": "[C;!R:2](=[O,S;!R])-[N;!R:1]-[#6]",
}
CREST_ENTROPY_OUTPUT_PARAMS: list[str] = [
    "S_conf",
    "S_rrho",
    "S_total",
]
THERMODYNAMIC_PARAMS: list[str] = [
    "E",  # electronic energy (incl. implicit solvation where used)
    "H",
    "S",
    "G",
]
ENERGY_MATRIX_PARAMS: list[str] = ["E", "H", "G"]

# CENSO sub-stage section names mapped to the file stem CENSO writes for them.
# Each section holds Boltzmann-averaged values; the matching ``<section>_min``
# section holds the values of the lowest free energy conformer.
CENSO_PARTS: dict[str, str] = {
    "censo_prescreening": "0_PRESCREENING",
    "censo_screening": "1_SCREENING",
    "censo_optimization": "2_OPTIMIZATION",
    "censo_refinement": "3_REFINEMENT",
}
CENSO_SECTIONS: list[str] = list(CENSO_PARTS)
MIN_SUFFIX = "_min"
PIPELINE_SECTIONS: list[str] = [
    "rdkit",
    "crest",
    *[
        f"{section}{suffix}"
        for section in CENSO_SECTIONS
        for suffix in ("", MIN_SUFFIX)
    ],
]
ENTHALPY_PLOTTING_SECTIONS: list[str] = [
    "rdkit",
    "crest",
    "censo_refinement",
]
ENTROPY_PLOTTING_SECTIONS: list[str] = [
    "rdkit",
    "crest",
    "censo_refinement",
    "censo_corrected",
]
