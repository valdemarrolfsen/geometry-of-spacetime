"""Independent algebraic checks; writes machine-readable results without sympy repr."""

import json
from pathlib import Path

import sympy as sp


def checks():
    eta, h, g, gp, lam = sp.symbols("eta h g gp lam", real=True)
    hp = -(h**2 + 1) / 2
    gpp = 2 * g - 3 * h * gp
    def diff(expr):
        return sp.diff(expr, h) * hp + sp.diff(expr, g) * gp + sp.diff(expr, gp) * gpp
    C = 2 / (3 * (h**2 + 1))
    V = -C * (gp + h * g)
    delta = C * (lam + 3) * g + 3 * h * V
    # The angular profile's Laplacian is represented by lam*F, without
    # restricting it to one actual eigenmode. Linearity allows superposition.
    identities = {
        "dust_euler": sp.simplify(diff(V) + h * V + g),
        "dust_continuity": sp.simplify(diff(delta) + lam * V - 3 * gp),
        "closed_dust_potential_coefficient": sp.simplify(2 * hp + h**2 - 1 + 2),
    }
    eps, F, Fchi, A = sp.symbols("epsilon F F_chi A", real=True)
    # A radial Euclidean graph X=A(1-eps F)n has induced radial coefficient
    # A^2[(1-eps F)^2+eps^2 F_chi^2]. Compare with the physical spatial metric.
    induced = A**2 * ((1 - eps * F)**2 + eps**2 * Fchi**2)
    identities["slice_embedding_first_order"] = sp.expand(induced - A**2 * (1 - 2 * eps * F)).coeff(eps, 1)
    v = sp.sqrt((1 + 2 * eps * F) / (1 - 2 * eps * F))
    identities["optical_speed_first_order"] = sp.simplify(sp.diff(v, eps).subs(eps, 0) - 2 * F)
    identities["lapse_first_order"] = sp.simplify(sp.diff(sp.sqrt(1 + 2 * eps * F), eps).subs(eps, 0) - F)
    result = {name: value == 0 for name, value in identities.items()}
    if not all(result.values()):
        raise AssertionError(identities)
    return result


if __name__ == "__main__":
    result = checks()
    root = Path(__file__).resolve().parents[1]
    (root / "results").mkdir(exist_ok=True)
    (root / "results/symbolic-checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{len(result)} symbolic identities verified.")
