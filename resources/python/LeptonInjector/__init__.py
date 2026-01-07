"""
LeptonInjector - Generates neutrino events for large volume Cherenkov telescopes.

This module provides tools for injecting neutrino events with proper cross-section
weighting for IceCube and similar detectors.

Basic usage:
    import LeptonInjector as LI

    # Create an injector
    injector = LI.Injector(
        NEvents=1000,
        FinalType1=LI.Particle.ParticleType.MuMinus,
        FinalType2=LI.Particle.ParticleType.Hadrons,
        DoublyDifferentialCrossSectionFile="path/to/xs.fits",
        TotalCrossSectionFile="path/to/xs_total.fits",
        Ranged=True
    )

    # Create controller and execute
    controller = LI.Controller(injector, ...)
    controller.Execute()
"""

import importlib as _importlib

# Import the compiled extension module
# The extension is named LeptonInjector.so and lives in the same directory
_ext_module = _importlib.import_module('.LeptonInjector', __name__)

# Export all public symbols from compiled extension
for _name in dir(_ext_module):
    if not _name.startswith('_'):
        globals()[_name] = getattr(_ext_module, _name)

# Clean up temporary variable
del _name

# Version - keep in sync with pyproject.toml and configure
__version__ = "1.1.0"

# Convenience: list commonly used classes for tab completion
__all__ = [
    # Core classes
    'Controller',
    'Injector',
    'Particle',
    'RNG',
    # Position/Direction
    'LI_Position',
    'LI_Direction',
    # Constants
    'Constants',
    # Helper functions
    'isLepton',
    'isCharged',
    'particleName',
    'particleMass',
    'kineticEnergy',
    'particleSpeed',
    'decideShape',
    'deduceInitialType',
    'getInteraction',
    # Coordinate functions
    'RotateX',
    'RotateY',
    'RotateZ',
    'rotateRelative',
    'computeCylinderIntersections',
]
