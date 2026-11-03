"""Analytical steady-state P1 transport in an infinite medium."""

from __future__ import annotations
import numpy as np
import MomentRT


class InfiniteP1:
    """Evaluate Phi = P exp(-sqrt(mu_a/D) r) / (4 pi D r).

    No boundary, Fresnel interface, or finite-source correction is applied.
    The observation window does not truncate the physical medium.

    Notes
    -----
    The Green function is Phi(r) = P exp(-sqrt(mu_a / D) r) / (4 pi D r).
    No physical outer boundary, Fresnel interface, or finite-source correction
    is applied. The observation window does not truncate the medium.
    """

    model_notes: tuple[str, ...] = (
        "InfiniteDiffusion is an analytical steady-state P1 Green function for an embedded isotropic point source in a homogeneous, unbounded medium.",
        "The point-source singularity is excluded; boundaries, finite source footprints, and refractive interfaces are not modeled.",
        "P1 diffusion is least reliable close to the source and where scattering does not dominate absorption.",
    )

    def solve(
        self,
        *,
        medium: MomentRT.medium.HomogeneousMedium,
        source: MomentRT.source.PointSource,
        mesh: MomentRT.mesh.CartesianMesh | MomentRT.mesh.TetrahedralMesh,
    ) -> MomentRT.fields.ScalarField:
        """Solve the transport equation for the supplied physical configuration.

        Parameters
        ----------
        medium : HomogeneousMedium
            Infinite homogeneous material resolved at the source wavelength.
        source : PointSource
            Embedded isotropic point source; no grid node may coincide with it.
        mesh : CartesianMesh
            Cartesian nodes used for sampling or finite-element discretization.

        Returns
        -------
        ScalarField
            Analytical infinite-medium fluence rate in watts per square meter.

        Raises
        ------
        TypeError
            If the medium or source is incompatible with the infinite model.
        ValueError
            If a grid node coincides with the singular point source.
        """

        if not isinstance(medium, MomentRT.medium.HomogeneousMedium):
            raise TypeError(
                "InfiniteDiffusion supports only HomogeneousMedium; finite structured samples require a boundary-aware solver."
            )

        if not isinstance(source, MomentRT.source.PointSource):
            raise TypeError(
                "InfiniteDiffusion requires an embedded PointSource, not surface illumination."
            )

        if not isinstance(mesh, MomentRT.mesh.CartesianMesh):
            raise TypeError("InfiniteDiffusion requires a Cartesian observation mesh.")

        medium = medium.at(wavelength=source.wavelength)

        absorption = medium.material.absorption.to("1/m").magnitude

        scattering = medium.material.reduced_scattering.to("1/m").magnitude

        diffusion_length = 1 / (3 * (absorption + scattering))

        displacement = mesh.points.to("m").magnitude - source.position.to("m").magnitude

        distance = np.linalg.norm(displacement, axis=-1)

        if np.any(distance == 0):
            raise ValueError(
                "The point source is singular; choose a grid that excludes its position."
            )

        values = (
            source.optical_power.to("W").magnitude
            * np.exp(-np.sqrt(absorption / diffusion_length) * distance)
            / (4 * np.pi * diffusion_length * distance)
        )

        return MomentRT.fields.ScalarField(
            mesh=mesh,
            values=values * MomentRT.units.ureg.watt / MomentRT.units.ureg.meter**2,
            name="Fluence rate",
            metadata={
                "transport": "infinite-medium steady-state P1 diffusion",
                "wavelength": source.wavelength,
                "source_position": source.position,
            },
        )
