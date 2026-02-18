#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Basic integration tests for the heat1d package."""

import pytest

import heat1d


def test_import():
    """Package imports successfully."""
    assert hasattr(heat1d, "Model")
    assert hasattr(heat1d, "Configurator")
    assert hasattr(heat1d, "Profile")


def test_version():
    """Version string is set."""
    assert heat1d.__version__ == "0.4.0"


def test_basic_model_run():
    """Basic model run completes without error."""
    m = heat1d.Model()
    m.run()
    assert m.T.shape[0] > 0
    assert m.T.shape[1] > 0
    # Surface temperature should be physically reasonable
    assert m.T[:, 0].max() > 300
    assert m.T[:, 0].min() > 50
