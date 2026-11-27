"""Analytic limits, independent quadrature, conservation laws, and source matching."""

import numpy as np
import pytest
from scipy.integrate import quad

from spacetime.model import ClosedDust, LocalizedPerturbation
from spacetime.rays import born_angular_shift, observer_frame, radial_emission_shift, shoot_to_source, trace_ray


@pytest.fixture(scope="module")
def model():
    return LocalizedPerturbation()


def test_closed_dust_friedmann_and_time():
    eta = np.linspace(0.1, 2.8, 100)
    b = ClosedDust()
    np.testing.assert_allclose(b.A(eta) * (b.h(eta)**2 + 1), 1, rtol=1e-13)
    step = 1e-5
    np.testing.assert_allclose((b.time(eta + step) - b.time(eta - step)) / (2 * step),
                               b.A(eta), atol=3e-11)


def test_growth_against_independent_integral(model):
    b = model.background
    def integral(eta):
        value = quad(lambda x: b.A(x) / b.h(x)**2, 0, eta, epsabs=1e-14)[0]
        return b.h(eta) * value / b.A(eta)**2
    normalization = integral(model.eta_observer)
    for eta in [0.1, 0.4, 0.9, 1.5, 1.8]:
        assert model.growth(eta) == pytest.approx(integral(eta) / normalization, rel=2e-9)


def test_global_compensation_and_low_modes(model):
    for power in [0, 1]:
        moment = quad(lambda chi: model.profile(np.cos(chi)) * np.cos(chi)**power
                      * np.sin(chi)**2, 0, np.pi, epsabs=1e-13)[0]
        assert abs(moment) < 1e-12
    integral = quad(lambda chi: model.density_comoving(1.8, np.cos(chi)) * np.sin(chi)**2,
                    0, np.pi, epsabs=1e-13)[0]
    assert abs(integral) < 1e-12
    assert model.profile(1) == pytest.approx(-1)
    assert model.density_comoving(1.8, 1) > 0
    assert min(model.density_comoving(1.8, np.cos(np.linspace(0, np.pi, 1000)))) < 0


def test_density_and_velocity_in_linear_regime(model):
    mu = np.cos(np.linspace(0, np.pi, 1001))
    for eta in [0.4, 1.0, 1.8]:
        assert np.max(np.abs(model.density_newtonian(eta, mu))) < 0.02
        assert np.max(np.abs(model.potential(eta, mu))) <= 1.001e-4


def test_profile_laplacian_independent_finite_difference(model):
    step = 2e-5
    for chi in [0.12, 0.3, 0.6, 1.3, 2.5]:
        f = lambda x: model.profile(np.cos(x))
        numerical = (f(chi + step) - 2*f(chi) + f(chi - step)) / step**2
        numerical += 2 / np.tan(chi) * (f(chi + step) - f(chi - step)) / (2*step)
        assert numerical == pytest.approx(model.laplacian_profile(np.cos(chi)), abs=2e-5)


def test_unperturbed_great_circle_and_redshift():
    model = LocalizedPerturbation(epsilon=0)
    ray = trace_ray(model, angle=0.3)
    n0, e0 = observer_frame(angle=0.3)
    expected = np.cos(ray["s"][:, None]) * n0 + np.sin(ray["s"][:, None]) * e0
    np.testing.assert_allclose(ray["n"], expected, atol=2e-9)
    expected_z = model.background.A(model.eta_observer) / model.background.A(ray["eta"]) - 1
    np.testing.assert_allclose(ray["redshift"], expected_z, atol=1e-8)


def test_ray_constraints_and_angular_momentum(model):
    ray = trace_ray(model, angle=0.25)
    n, p = ray["n"], ray["p"]
    np.testing.assert_allclose(np.sum(n*n, axis=1), 1, atol=2e-9)
    np.testing.assert_allclose(np.sum(n*p, axis=1), 0, atol=2e-9)
    angular = n[:, 1]*p[:, 2] - n[:, 2]*p[:, 1]
    np.testing.assert_allclose(angular, angular[0], atol=2e-9)
    # Test the null condition with an independent finite difference of the path.
    dn = (n[2:] - n[:-2]) / (ray["s"][2:, None] - ray["s"][:-2, None])
    phi = model.potential(ray["eta"][1:-1], n[1:-1, 0])
    null_residual = -(1+2*phi) + (1-2*phi)*np.sum(dn*dn, axis=1)
    assert np.max(np.abs(null_residual)) < 3e-6


def test_born_delay_and_second_order_remainder():
    errors = []
    for eps in [1e-4, 2e-4]:
        result = radial_emission_shift(LocalizedPerturbation(epsilon=eps))
        assert result["additional_conformal_lookback"] > 0
        assert result["additional_conformal_lookback"] == pytest.approx(result["born_conformal_lookback"], rel=1e-3)
        errors.append(abs(result["additional_conformal_lookback"] - result["born_conformal_lookback"]))
    assert 3.7 < errors[1] / errors[0] < 4.3


def test_matched_source_lensing_and_reflection(model):
    zero = shoot_to_source(LocalizedPerturbation(epsilon=0))
    assert abs(zero["angular_shift_rad"]) < 1e-8
    positive = shoot_to_source(model, baseline_angle=0.25)
    negative = shoot_to_source(model, baseline_angle=-0.25)
    assert positive["source_residual"] < 2e-8
    assert positive["angular_shift_rad"] > 0
    assert positive["angular_shift_rad"] == pytest.approx(-negative["angular_shift_rad"], abs=2e-9)


def test_solver_tolerance_stability(model):
    normal = trace_ray(model, angle=0.25)
    tight = trace_ray(model, angle=0.25, rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(normal["n"][-1], tight["n"][-1], atol=2e-9)
    assert normal["redshift"][-1] == pytest.approx(tight["redshift"][-1], abs=2e-8)


def test_lensing_against_independent_first_order_integral(model):
    result = shoot_to_source(model)
    assert result["angular_shift_rad"] == pytest.approx(born_angular_shift(model), rel=0.002)


def test_frequency_against_endpoint_and_integrated_potential(model):
    ray = trace_ray(model, angle=0.25, samples=4001, rtol=2e-12, atol=2e-14)
    phi = model.potential(ray["eta"], ray["n"][:, 0])
    phi_prime = model.potential_prime(ray["eta"], ray["n"][:, 0])
    integrated = np.trapezoid(2*phi_prime/(1-4*phi**2), ray["s"])
    expected = (model.background.A(model.eta_observer) / model.background.A(ray["eta"][-1])
                * np.sqrt((1+2*phi[0])/(1+2*phi[-1])) * np.exp(-integrated))
    assert ray["redshift"][-1]+1 == pytest.approx(expected, rel=3e-10)


def test_invalid_intervals_and_model_parameters():
    with pytest.raises(ValueError):
        LocalizedPerturbation(eta_observer=3.5)
    with pytest.raises(ValueError):
        LocalizedPerturbation(epsilon=0.1)
    with pytest.raises(ValueError):
        trace_ray(LocalizedPerturbation(), lookback=2)
