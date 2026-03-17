"""
Unicode math formatting utilities for the Alpha Ladder dashboard.

Provides Greek letters, superscripts, subscripts, and symbol helpers
for replacing ASCII math notation with proper Unicode in display strings.
"""

GREEK = {
    "alpha": "\u03b1", "beta": "\u03b2", "gamma": "\u03b3", "delta": "\u03b4",
    "epsilon": "\u03b5", "zeta": "\u03b6", "eta": "\u03b7", "theta": "\u03b8",
    "kappa": "\u03ba", "lambda": "\u03bb", "mu": "\u03bc", "nu": "\u03bd",
    "pi": "\u03c0", "rho": "\u03c1", "sigma": "\u03c3", "tau": "\u03c4",
    "phi": "\u03c6", "chi": "\u03c7", "psi": "\u03c8", "omega": "\u03c9",
    "Lambda": "\u039b", "Sigma": "\u03a3", "Omega": "\u03a9",
}

SUP = {
    "0": "\u2070", "1": "\u00b9", "2": "\u00b2", "3": "\u00b3",
    "4": "\u2074", "5": "\u2075", "6": "\u2076", "7": "\u2077",
    "8": "\u2078", "9": "\u2079", "-": "\u207b", "+": "\u207a",
}

SUB = {
    "0": "\u2080", "1": "\u2081", "2": "\u2082", "3": "\u2083",
    "4": "\u2084", "5": "\u2085", "6": "\u2086", "7": "\u2087",
    "8": "\u2088", "9": "\u2089",
}

SYMBOLS = {
    "hbar": "\u210f", "sqrt": "\u221a", "approx": "\u2248", "pm": "\u00b1",
    "times": "\u00d7", "cdot": "\u00b7", "arrow": "\u2192", "leq": "\u2264",
    "geq": "\u2265", "neq": "\u2260", "inf": "\u221e", "partial": "\u2202",
    "nabla": "\u2207", "ell": "\u2113",
}


def sup(n):
    """Convert a string or int to Unicode superscript. sup(24) -> '\u00b2\u2074'."""
    return "".join(SUP.get(c, c) for c in str(n))


def sub(n):
    """Convert a string or int to Unicode subscript. sub(0) -> '\u2080'."""
    return "".join(SUB.get(c, c) for c in str(n))
