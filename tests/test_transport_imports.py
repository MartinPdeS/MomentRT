"""The transport package exposes the solver classes at the package root."""

import MomentRT
import MomentRT.transport.diffusion
import MomentRT.transport.p3
import MomentRT.transport.sp3


def test_solver_import_paths_share_classes():
    assert MomentRT.Diffusion is MomentRT.transport.diffusion.Diffusion
    assert MomentRT.InfiniteDiffusion is MomentRT.transport.diffusion.InfiniteDiffusion
    assert MomentRT.SP3 is MomentRT.transport.sp3.SP3
    assert MomentRT.P3 is MomentRT.transport.p3.P3
