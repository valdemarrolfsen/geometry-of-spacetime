"""Past-directed null rays on an equatorial S^2 of the three-sphere.

Integrates geodesics of the specified truncated metric, exactly in epsilon at
the numerical level. Physical accuracy remains first order in epsilon.
The ray parameter s=eta_o-eta is conformal lookback, not affine parameter.
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import root


def observer_frame(chi=0.65, angle=0.0):
    n = np.array([np.cos(chi), np.sin(chi), 0.0])
    toward_lens = np.array([np.sin(chi), -np.cos(chi), 0.0])
    e = np.cos(angle) * toward_lens + np.sin(angle) * np.array([0.0, 0.0, 1.0])
    return n, e


def trace_ray(model, *, chi_observer=0.65, angle=0.0, lookback=1.3,
              samples=501, rtol=2e-10, atol=2e-12):
    if not 0 < lookback < model.eta_observer - model.eta_initial:
        raise ValueError("Ray must remain inside the background's expanding integration interval.")
    n0, e0 = observer_frame(chi_observer, angle)

    def rhs(s, y):
        n, p = y[:3], y[3:]
        pnorm = np.linalg.norm(p)
        eta = model.eta_observer - s
        v = model.optical_speed(eta, n[0])
        gradv = model.optical_gradient(eta, n)
        return np.concatenate([v * p / pnorm, -pnorm * gradv - v * pnorm * n])

    sol = solve_ivp(rhs, (0, lookback), np.concatenate([n0, e0]),
                    t_eval=np.linspace(0, lookback, samples), method="DOP853",
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    n, p = sol.y[:3].T, sol.y[3:].T
    eta = model.eta_observer - sol.t
    # Frequency measured by the normal observers; dust Doppler terms excluded.
    energy = np.linalg.norm(p, axis=1) / (model.background.A(eta) * model.spatial_factor(eta, n[:, 0]))
    return {"s": sol.t, "eta": eta, "n": n, "p": p,
            "redshift": energy / energy[0] - 1}


def radial_emission_shift(model, *, chi_observer=0.65, distance=1.3):
    """Fixed reception eta and fixed comoving emitter on the opposite side.

Returns additional conformal lookback and positive earlier-emission cosmic
time shift in units A_max/c, against the matched FLRW reference. It is a
gauge-specified comparison, not a raw coordinate-independent observable.
"""
    max_s = model.eta_observer - model.eta_initial
    if not 0 < distance < max_s:
        raise ValueError("Source must fit inside the modeled conformal-time interval.")

    def rhs(s, y):
        return [model.optical_speed(model.eta_observer - s, np.cos(chi_observer - y[0]))]

    def event(s, y):
        return y[0] - distance

    event.terminal, event.direction = True, 1
    sol = solve_ivp(rhs, (0, max_s), [0.0], events=event, method="DOP853",
                    rtol=2e-12, atol=2e-14, max_step=0.03)
    if not sol.success or len(sol.t_events[0]) != 1:
        raise RuntimeError("The source was not reached in the modeled interval.")
    s = float(sol.t_events[0][0])
    born = -2 * quad(lambda ell: model.potential(model.eta_observer - ell,
                                                np.cos(chi_observer - ell)),
                     0, distance, epsabs=1e-14)[0]
    eta_background = model.eta_observer - distance
    return {"additional_conformal_lookback": s - distance,
            "born_conformal_lookback": born,
            "earlier_emission_time_Amax_over_c": float(model.background.time(eta_background)
                                                     - model.background.time(model.eta_observer - s)),
            "emission_eta": model.eta_observer - s}


def shoot_to_source(model, *, baseline_angle=0.25, baseline_distance=1.3,
                    chi_observer=0.65):
    """Match a coordinate-stationary source worldline and reception event.

    The change in launch angle is the change in the receiver's orthonormal
    frame, since the spatial metric is conformally isotropic there. This is
    a model-to-baseline astrometric prediction with an explicit source match.
    """
    n0, e0 = observer_frame(chi_observer, baseline_angle)
    source = n0 * np.cos(baseline_distance) + e0 * np.sin(baseline_distance)
    basis1 = -n0 * np.sin(baseline_distance) + e0 * np.cos(baseline_distance)
    basis2 = np.cross(source, basis1)

    def residual(values):
        angle, distance = values
        if not 0 < distance < model.eta_observer - model.eta_initial:
            raise ValueError("Shooting left the modeled expanding interval.")
        ray = trace_ray(model, chi_observer=chi_observer, angle=angle,
                        lookback=distance, samples=2, rtol=2e-11, atol=2e-13)
        difference = ray["n"][-1] - source
        return [np.dot(difference, basis1), np.dot(difference, basis2)]

    solution = root(residual, [baseline_angle, baseline_distance], tol=1e-9)
    error = np.linalg.norm(residual(solution.x))
    if error > 2e-9:
        raise RuntimeError(f"Source shooting did not converge: {solution.message}; residual={error}")
    angle, distance = map(float, solution.x)
    ray = trace_ray(model, chi_observer=chi_observer, angle=angle,
                    lookback=distance)
    return {"angle": angle, "angular_shift_rad": angle - baseline_angle,
            "additional_conformal_lookback": distance - baseline_distance,
            "source_residual": float(np.linalg.norm(ray["n"][-1] - source)),
            "source": source, "ray": ray}


def born_angular_shift(model, *, baseline_angle=0.25, baseline_distance=1.3,
                       chi_observer=0.65):
    """First-order fixed-source angle from the sphere's transverse Jacobi equation."""
    n0, e0 = observer_frame(chi_observer, baseline_angle)
    _, perpendicular = observer_frame(chi_observer, baseline_angle + np.pi / 2)

    def integrand(s):
        n = n0*np.cos(s) + e0*np.sin(s)
        grad = model.potential_gradient(model.eta_observer-s, n)
        return np.sin(baseline_distance-s)*np.dot(grad, perpendicular)

    return 2*quad(integrand, 0, baseline_distance, epsabs=1e-14)[0] / np.sin(baseline_distance)
