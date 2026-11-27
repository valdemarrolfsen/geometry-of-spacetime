import numpy as np

from spacetime.diagrams import event_coordinates, inverse_nested_expanding, visible_ring
from spacetime.model import ClosedDust


def test_maps_recover_the_same_events_on_expanding_branch():
    eta = np.linspace(0.05, 3.0, 21)
    chi = np.linspace(-2.8, 2.8, 21)
    events = event_coordinates(eta, chi)
    recovered_eta, recovered_chi = inverse_nested_expanding(*events["nested"].T)
    np.testing.assert_allclose(recovered_eta, eta, atol=2e-14)
    np.testing.assert_allclose(recovered_chi, chi, atol=2e-14)
    np.testing.assert_allclose(events["proper"][:, 0],
                               ClosedDust.A(eta)*events["conformal"][:, 0])


def test_recollapse_overlaps_distinct_events_in_nested_map():
    events = event_coordinates([np.pi/2, 3*np.pi/2], [0.8, 0.8])
    np.testing.assert_allclose(events["nested"][0], events["nested"][1], atol=1e-14)
    assert events["comoving"][1, 1] > events["comoving"][0, 1] + 2


def test_observed_ring_is_a_subset_of_a_spatial_slice():
    eta, eo = 0.9, 1.8
    ring = visible_ring(eta, np.linspace(0, 2*np.pi, 101), eo)
    radius = ClosedDust.A(eta)/ClosedDust.A(eo)
    np.testing.assert_allclose(np.linalg.norm(ring, axis=1), radius)
    separation = np.arccos(ring[:, 0]/radius)
    np.testing.assert_allclose(separation, eo-eta)
