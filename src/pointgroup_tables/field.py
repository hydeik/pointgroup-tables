"""The document field ``Q(zeta_M)`` (docs/spec/values.md, docs/spec/basis-conventions.md).

Only the metadata is needed here: the canonical conductor of a group and the
degree ``phi(M)``. The arithmetic itself belongs to the values chapter.
"""

from __future__ import annotations

import math
from typing import Any

from .rotation import RotationGroup

POLYHEDRAL_CONDUCTOR = {"T": 12, "O": 8, "I": 20}


def euler_phi(m: int) -> int:
    return sum(1 for k in range(1, m + 1) if math.gcd(k, m) == 1)


def canonical_conductor(R: RotationGroup) -> int:
    """The canonical conductor, fixed by ``R`` alone."""
    if R.axis_order is None:
        return POLYHEDRAL_CONDUCTOR[R.name]
    N = R.axis_order
    if R.name.startswith("D"):
        return math.lcm(2 * N, 4)
    return N if N % 2 else 2 * N


def field_metadata(conductor: int) -> dict[str, Any]:
    return {
        "kind": "cyclotomic",
        "conductor": conductor,
        "degree": euler_phi(conductor),
        "generator": f"exp(2*pi*i/{conductor})",
        "basis": "power",
    }
