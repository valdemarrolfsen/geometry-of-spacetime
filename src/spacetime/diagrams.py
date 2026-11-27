"""Coordinate maps for a controlled comparison of cosmological diagrams.

All maps represent the same closed-dust events. The nested map is a drawing
of spatial slices, not a Euclidean isometry of Lorentzian spacetime.
The full dust history 0 < eta < 2*pi is allowed to expose its time folding.
"""

import numpy as np

from spacetime.model import ClosedDust


def event_coordinates(eta, chi, eta_observer=1.8):
    eta, chi = np.broadcast_arrays(np.asarray(eta), np.asarray(chi))
    if np.any((eta <= 0) | (eta >= 2*np.pi)):
        raise ValueError("Use the regular closed-dust interval 0 < eta < 2*pi.")
    if not 0 < eta_observer < 2*np.pi:
        raise ValueError("Reception must be inside the regular dust history.")
    b = ClosedDust()
    radius = b.A(eta) / b.A(eta_observer)
    time = b.time(eta)
    return {
        "comoving": np.stack([chi, time], axis=-1),
        "proper": np.stack([b.A(eta)*chi, time], axis=-1),
        "conformal": np.stack([chi, eta], axis=-1),
        "nested": np.stack([radius*np.cos(chi), radius*np.sin(chi)], axis=-1),
    }


def inverse_nested_expanding(x, y, eta_observer=1.8):
    """Invert only the expanding branch; chi is returned in [-pi, pi]."""
    radius = np.hypot(x, y)
    scale_factor = radius * ClosedDust.A(eta_observer)
    if np.any((scale_factor <= 0) | (scale_factor >= 1)):
        raise ValueError("The inverse requires 0 < A/A_max < 1.")
    return 2*np.arcsin(np.sqrt(scale_factor)), np.arctan2(y, x)


def visible_ring(eta, angle, eta_observer=1.8):
    """Past-cone intersection with a totally geodesic S^2 spatial section.

The full spatial slice is S^3. Its full past-cone intersection is S^2;
the circle returned here is one section of that two-sphere.
"""
    chi = eta_observer - eta
    if not 0 < chi < np.pi:
        raise ValueError("Use the direct past cone before the antipodal caustic.")
    angle = np.asarray(angle)
    radius = ClosedDust.A(eta) / ClosedDust.A(eta_observer)
    return radius * np.stack([np.full_like(angle, np.cos(chi)),
                              np.sin(chi)*np.cos(angle),
                              np.sin(chi)*np.sin(angle)], axis=-1)
