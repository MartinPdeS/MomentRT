"""The transport package exposes the solver classes at the package root."""

import MomentRT
import MomentRT.transport.diffusion
import MomentRT.transport.infinite_p1
import MomentRT.transport.p1
import MomentRT.transport.p3
import MomentRT.transport.sp3
import MomentRT.transport.transient_p1


def test_solver_import_paths_share_classes():
    assert MomentRT.P1 is MomentRT.transport.p1.P1
    assert MomentRT.InfiniteP1 is MomentRT.transport.infinite_p1.InfiniteP1
    assert MomentRT.Diffusion is MomentRT.P1
    assert MomentRT.InfiniteDiffusion is MomentRT.InfiniteP1
    assert MomentRT.Diffusion is MomentRT.transport.diffusion.Diffusion
    assert MomentRT.InfiniteDiffusion is MomentRT.transport.diffusion.InfiniteDiffusion
    assert MomentRT.SP3 is MomentRT.transport.sp3.SP3
    assert MomentRT.P3 is MomentRT.transport.p3.P3
    assert MomentRT.transport.P1 is MomentRT.P1
    assert MomentRT.transport.InfiniteP1 is MomentRT.InfiniteP1
    assert MomentRT.transport.SP3 is MomentRT.SP3
    assert MomentRT.transport.P3 is MomentRT.P3
    assert MomentRT.transport.TransientP1 is MomentRT.TransientP1
    assert MomentRT.transport.Diffusion is MomentRT.Diffusion
    assert MomentRT.transport.InfiniteDiffusion is MomentRT.InfiniteDiffusion
