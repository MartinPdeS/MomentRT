"""Ideal isotropic point illumination inside an infinite medium."""

from .units import validate_units


class PointSource:
    def __init__(self, *, position, wavelength, optical_power):
        validate_units(value=position, unit="m", name="position", shape=(3,))

        validate_units(
            value=wavelength, unit="m", name="wavelength", shape=(), positive=True
        )

        validate_units(
            value=optical_power, unit="W", name="optical_power", shape=(), positive=True
        )

        self.position = position

        self.wavelength = wavelength

        self.optical_power = optical_power


class UniformSurfaceSource:
    """Uniform incoming diffuse partial current across one complete sample face.

    optical_power is power already transmitted into the sample, before diffuse
    internal reflection. This does not model collimated LED entry or specular
    Fresnel loss. The boundary condition's source factor is applied by the solver.
    """

    def __init__(self, *, face, wavelength, optical_power):
        from .boundary import FACES

        if face not in FACES:
            raise ValueError(f"face must be one of {FACES}.")

        validate_units(
            value=wavelength, unit="m", name="wavelength", shape=(), positive=True
        )

        validate_units(
            value=optical_power, unit="W", name="optical_power", shape=(), positive=True
        )

        self.face = face

        self.wavelength = wavelength

        self.optical_power = optical_power


class SurfaceSource(UniformSurfaceSource):
    """Circular top-hat incoming diffuse illumination on an axis-aligned face.

    position is an absolute 3D surface coordinate; radius is the spot radius.
    The entire disk must fit on the selected face. Surface quadrature samples
    its footprint and is normalized to optical_power; refine the mesh to
    resolve the disk edge and check the reported quadrature spot area.
    """

    def __init__(self, *, face, position, radius, wavelength, optical_power):
        super().__init__(face=face, wavelength=wavelength, optical_power=optical_power)

        validate_units(value=position, unit="m", name="position", shape=(3,))

        validate_units(value=radius, unit="m", name="radius", shape=(), positive=True)

        self.position = position.copy()

        self.radius = radius.copy()
