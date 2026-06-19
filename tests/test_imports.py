"""Verify that all package modules import cleanly."""


def test_import_package():
    import cosmic_string_signals
    assert hasattr(cosmic_string_signals, "__version__")


def test_import_geometry():
    from cosmic_string_signals import geometry
    assert hasattr(geometry, "rotation_x")
    assert hasattr(geometry, "rotation_y")
    assert hasattr(geometry, "rotation_z")
    assert hasattr(geometry, "apply_rotation")
    assert hasattr(geometry, "build_wake_wedge")
    assert hasattr(geometry, "point_in_hull")


def test_import_temperature():
    from cosmic_string_signals import temperature
    assert hasattr(temperature, "brightness_temperature")
    assert hasattr(temperature, "photon_cmb_temperature")
    assert hasattr(temperature, "kinetic_temperature")
    assert hasattr(temperature, "gas_temperature")


def test_import_filtering():
    from cosmic_string_signals import filtering
    assert hasattr(filtering, "unfold")
    assert hasattr(filtering, "matched_filter_1d")
    assert hasattr(filtering, "matched_filter_2d")
