"""Classical single-phase pipe correlations used as validation references."""
from __future__ import annotations

import math


def petukhov_friction(re: float) -> float:
    """Darcy friction factor for a smooth pipe (Petukhov), 3e3 < Re < 5e6."""
    return (0.79 * math.log(re) - 1.64) ** -2


def gnielinski(re: float, pr: float) -> float:
    """Nusselt number, Gnielinski (1976), 3e3 < Re < 5e6, 0.5 < Pr < 2000."""
    f = petukhov_friction(re)
    return (f / 8.0) * (re - 1000.0) * pr / (
        1.0 + 12.7 * math.sqrt(f / 8.0) * (pr ** (2.0 / 3.0) - 1.0)
    )


def dittus_boelter(re: float, pr: float, heating: bool = True) -> float:
    """Nusselt number, Dittus-Boelter: 0.023 Re^0.8 Pr^n, n = 0.4 (heating) or 0.3."""
    n = 0.4 if heating else 0.3
    return 0.023 * re ** 0.8 * pr ** n
