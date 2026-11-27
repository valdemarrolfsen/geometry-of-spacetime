"""A smooth compensated scalar perturbation of closed dust FLRW.

The source is reconstructed from the linear Einstein constraints. It is an
evolving dust overdensity, not a fixed Schwarzschild mass or a black hole.
Spatial coordinates are unit vectors on S^3 (or a totally geodesic S^2).
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ive


@dataclass(frozen=True)
class ClosedDust:
    """Expanding branch, 0 < eta < pi; A and ct in units of A_max."""

    @staticmethod
    def A(eta):
        return np.sin(np.asarray(eta) / 2) ** 2

    @staticmethod
    def time(eta):
        eta = np.asarray(eta)
        return (eta - np.sin(eta)) / 2

    @staticmethod
    def h(eta):
        return 1 / np.tan(np.asarray(eta) / 2)

    @staticmethod
    def h_prime(eta):
        return -0.5 / np.sin(np.asarray(eta) / 2) ** 2


@dataclass
class LocalizedPerturbation:
    epsilon: float = 1e-4
    sigma: float = 0.18
    eta_observer: float = 1.8
    eta_initial: float = 1e-3
    background: ClosedDust = field(default_factory=ClosedDust)

    def __post_init__(self):
        if not 0 < self.eta_initial < self.eta_observer < np.pi:
            raise ValueError("Use a strictly expanding interval 0 < eta_i < eta_o < pi.")
        if not 0 < self.sigma <= 0.5 or abs(self.epsilon) >= 0.01:
            raise ValueError("Require 0 < sigma <= 0.5 and |epsilon| < 0.01; check density separately.")
        self.kappa = self.sigma ** -2
        self.q0 = 2 * ive(1, self.kappa) / self.kappa
        self.q1 = 8 * ive(2, self.kappa) / self.kappa
        self.normalization = 1 - self.q0 - self.q1
        e = self.eta_initial
        initial = [1 + e**2 / 7 + e**4 / 84 + e**6 / 1320,
                   2 * e / 7 + e**3 / 21 + e**5 / 220]
        self._growth = solve_ivp(
            lambda eta, y: [y[1], 2 * y[0] - 3 * self.background.h(eta) * y[1]],
            (e, self.eta_observer), initial, dense_output=True,
            method="DOP853", rtol=2e-12, atol=2e-14,
        )
        if not self._growth.success:
            raise RuntimeError(self._growth.message)
        self._growth_norm = self._growth.sol(self.eta_observer)[0]

    def growth(self, eta, derivative=False):
        eta = np.asarray(eta)
        if np.any(eta < self.eta_initial) or np.any(eta > self.eta_observer + 1e-12):
            raise ValueError("Growth requested outside the integrated conformal-time interval.")
        return self._growth.sol(eta)[int(derivative)] / self._growth_norm

    def profile(self, mu, derivative=0):
        """F(cos chi), with S^3 monopole and dipole removed and F(1)=-1."""
        mu = np.asarray(mu)
        q = np.exp(self.kappa * (mu - 1))
        if derivative == 0:
            return -(q - self.q0 - self.q1 * mu) / self.normalization
        if derivative == 1:
            return -(self.kappa * q - self.q1) / self.normalization
        if derivative == 2:
            return -self.kappa**2 * q / self.normalization
        raise ValueError("Only derivative orders 0, 1, 2 are available.")

    def laplacian_profile(self, mu):
        mu = np.asarray(mu)
        return (1 - mu**2) * self.profile(mu, 2) - 3 * mu * self.profile(mu, 1)

    def potential(self, eta, mu):
        return self.epsilon * self.growth(eta) * self.profile(mu)

    def potential_prime(self, eta, mu):
        return self.epsilon * self.growth(eta, True) * self.profile(mu)

    def potential_gradient(self, eta, n):
        """Tangential S^d gradient; the lens is at (1,0,...,0)."""
        n = np.asarray(n)
        lens = np.zeros_like(n)
        lens[0] = 1
        return self.epsilon * self.growth(eta) * self.profile(n[0], 1) * (lens - n[0] * n)

    def density_comoving(self, eta, mu):
        h = self.background.h(eta)
        return (2 * self.epsilon * self.growth(eta) / (3 * (h**2 + 1))
                * (self.laplacian_profile(mu) + 3 * self.profile(mu)))

    def velocity_potential(self, eta, mu):
        h = self.background.h(eta)
        return -2 * (self.potential_prime(eta, mu) + h * self.potential(eta, mu)) / (3 * (h**2 + 1))

    def density_newtonian(self, eta, mu):
        return self.density_comoving(eta, mu) + 3 * self.background.h(eta) * self.velocity_potential(eta, mu)

    def lapse(self, eta, mu):
        """Exact square root of the chosen truncated metric, not an exact GR solution."""
        return np.sqrt(1 + 2 * self.potential(eta, mu))

    def spatial_factor(self, eta, mu):
        return np.sqrt(1 - 2 * self.potential(eta, mu))

    def optical_speed(self, eta, mu):
        return self.lapse(eta, mu) / self.spatial_factor(eta, mu)

    def optical_gradient(self, eta, n):
        phi = self.potential(eta, n[0])
        return (2 * self.optical_speed(eta, n[0]) / (1 - 4 * phi**2)
                * self.potential_gradient(eta, n))

    def curvature_perturbation(self, eta, mu):
        """A^2 delta(3-R), a dimensionless first-order curvature perturbation."""
        return 4 * self.epsilon * self.growth(eta) * (self.laplacian_profile(mu) + 3 * self.profile(mu))
