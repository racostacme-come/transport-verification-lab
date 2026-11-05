"""Independent constitutive, event-boundary and energy accounting checks."""
from dataclasses import FrozenInstanceError, replace
import json
import numpy as np
import pytest
from phasechange import Material, integrate_uniform
from phasechange.reference import event_temperature
from phasechange.__main__ import ROOT, experiment, check


@pytest.fixture
def material():
    return Material(2., 3., 5., 12., 10., 14., 4.)


def test_hand_calculated_enthalpy_and_inverse(material):
    t = [0., 4., 10., 11., 12., 14., 20.]
    h = [-24., 0., 36., 48., 60., 84., 120.]
    np.testing.assert_array_equal(material.enthalpy(t), h)
    np.testing.assert_array_equal(material.temperature(h), t)
    np.testing.assert_array_equal(material.liquid_fraction(t), [0, 0, 0, .25, .5, 1, 1])


def test_monotone_and_two_way_roundtrip(material):
    # Include adjacent floating point numbers around both phase kinks.
    t = np.unique(np.r_[np.linspace(0, 100, 1001), 10., 14.,
                        np.nextafter([10., 14.], -np.inf), np.nextafter([10., 14.], np.inf)])
    h = material.enthalpy(t)
    assert np.all(np.diff(h) > 0)
    np.testing.assert_allclose(material.temperature(h), t, rtol=2e-15, atol=2e-14)
    trial = np.linspace(-24, 1000, 1001)
    np.testing.assert_allclose(material.enthalpy(material.temperature(trial)), trial, rtol=2e-14, atol=2e-13)


@pytest.mark.parametrize('temperature', [8., 11., 16.])
def test_centered_derivative_away_from_kinks(material, temperature):
    eps = 1e-4
    numerical = (material.enthalpy(temperature+eps)-material.enthalpy(temperature-eps))/(2*eps)
    assert numerical == pytest.approx(float(material.tangent(temperature)), rel=2e-10)


@pytest.mark.parametrize('temperature,left,right', [(10., 6., 12.), (14., 12., 6.)])
def test_one_sided_endpoint_derivatives(material, temperature, left, right):
    eps = 1e-4
    assert material.tangent(temperature, side='left') == left
    assert material.tangent(temperature, side='right') == right
    assert (material.enthalpy(temperature)-material.enthalpy(temperature-eps))/eps == pytest.approx(left)
    assert (material.enthalpy(temperature+eps)-material.enthalpy(temperature))/eps == pytest.approx(right)


def test_zero_latent_limit_is_linear_everywhere(material):
    m = replace(material, latent_heat=0.)
    t = np.linspace(0, 20, 201)
    np.testing.assert_allclose(m.enthalpy(t), 6*(t-4))
    np.testing.assert_allclose(m.temperature(6*(t-4)), t, atol=2e-15)
    np.testing.assert_array_equal(m.tangent(t), np.full_like(t, 6.))


def test_reference_offset_changes_no_temperature_or_energy_increment(material):
    shifted = replace(material, reference=100.)
    a = integrate_uniform(material, 8., [3., 5., 2.], [12., 12., -12.])
    b = integrate_uniform(shifted, 8., [3., 5., 2.], [12., 12., -12.])
    np.testing.assert_allclose(a.temperature, b.temperature)
    np.testing.assert_allclose(np.asarray(b.enthalpy)-a.enthalpy, -576.)
    np.testing.assert_allclose(b.energy_defect, 0., atol=1e-12)


@pytest.mark.parametrize('field,value', [('density', 0), ('heat_capacity', -1), ('conductivity', 0),
    ('latent_heat', -1), ('solidus', -1), ('liquidus', 10), ('reference', -1),
    ('density', float('nan')), ('latent_heat', float('inf')), ('density', 1e308)])
def test_invalid_material_rejected(material, field, value):
    with pytest.raises(ValueError):
        replace(material, **{field: value})


@pytest.mark.parametrize('method,value', [('enthalpy', -1), ('enthalpy', np.nan),
    ('temperature', -25), ('temperature', np.inf), ('liquid_fraction', [-1, 12]),
    ('tangent', np.nan), ('enthalpy', 1e308)])
def test_invalid_state_rejected(material, method, value):
    with pytest.raises(ValueError):
        getattr(material, method)(value)


def test_immutable_contracts(material):
    with pytest.raises(FrozenInstanceError):
        material.density = 4
    data = integrate_uniform(material, 8, [1], [12])
    with pytest.raises(TypeError):
        data.temperature[0] = 3


@pytest.mark.parametrize('durations,powers', [([], []), ([0], [1]), ([-1], [1]),
    ([1, 2], [1]), ([[1]], [[1]]), ([np.inf], [1]), ([1], [np.nan]), ([1e308], [1e308])])
def test_invalid_heating_schedule_rejected(material, durations, powers):
    with pytest.raises(ValueError):
        integrate_uniform(material, 8., durations, powers)


def test_uniform_crosses_both_endpoints_and_reverses_exactly(material):
    result = integrate_uniform(material, 8., [10., 10.], [12., -12.])
    np.testing.assert_array_equal(result.temperature, [8., 24., 8.])
    np.testing.assert_array_equal(result.supplied_heat, [0., 120., 0.])
    np.testing.assert_array_equal(result.energy_defect, [0., 0., 0.])


@pytest.mark.parametrize('initial,duration,power,expected', [
    (8, 1, 12, 10), (10, 4, 12, 14), (14, 1, -12, 13),
    (14, 4, -12, 10), (10, 1, -12, 8), (8, 10, 12, 24), (24, 10, -12, 8),
    (12, 10, 0, 12), (12, 0, 12, 12)])
def test_independent_event_reference(initial, duration, power, expected):
    actual = event_temperature(initial, duration, power, density=2., heat_capacity=3.,
                               latent_heat=12., solidus=10., liquidus=14.)
    assert actual == expected


def test_unphysical_cooling_is_not_clipped(material):
    with pytest.raises(ValueError):
        integrate_uniform(material, 8, [10], [-120])
    with pytest.raises(ValueError):
        event_temperature(8, 10, -120, density=2, heat_capacity=3, latent_heat=12, solidus=10, liquidus=14)


def test_fixed_experiment_and_versioned_evidence():
    result = experiment(); check(result)
    recorded = json.loads((ROOT/'results/enthalpy-baseline.json').read_text(encoding='utf-8'))
    assert result['source_sha256'] == recorded['source_sha256']
    assert len(result['cases']) == 30
    for a, b in zip(result['cases'], recorded['cases']):
        assert (a['case'], a['subdivisions']) == (b['case'], b['subdivisions'])
        np.testing.assert_allclose(a['history']['temperature'], b['history']['temperature'], rtol=0, atol=1e-11)
        assert a['scaled_energy_defect'] <= 1e-12
    cycles = [x for x in result['cases'] if x['case']=='heat_cool_cycle']
    assert all(abs(x['history']['temperature'][-1]-305.) < 1e-11 for x in cycles)


def test_acceptance_rejects_nonfinite_and_large_error():
    result = experiment()
    for value in (np.nan, 1e-6):
        with pytest.raises(RuntimeError):
            check({**result, 'max_scaled_temperature_error': value})
