"""Phase-function modules preserve the existing public class identities."""

import pickle
import MomentRT
from MomentRT.phase_function.base import PhaseFunction
from MomentRT.phase_function.henyey_greenstein import HenyeyGreenstein
from MomentRT.phase_function.isotropic import Isotropic
from MomentRT.phase_function.rayleigh import Rayleigh
from MomentRT.phase_function.tabulated import Tabulated


def test_phase_function_import_paths_share_classes() -> None:
    """Each direct module and convenience export names the same class."""
    assert (
        MomentRT.PhaseFunction is MomentRT.phase_function.PhaseFunction is PhaseFunction
    )
    assert (
        MomentRT.HenyeyGreenstein
        is MomentRT.phase_function.HenyeyGreenstein
        is HenyeyGreenstein
    )
    assert MomentRT.Isotropic is MomentRT.phase_function.Isotropic is Isotropic
    assert MomentRT.Rayleigh is MomentRT.phase_function.Rayleigh is Rayleigh
    assert (
        MomentRT.TabulatedPhaseFunction
        is MomentRT.phase_function.Tabulated
        is Tabulated
    )


def test_old_pickle_class_path_still_resolves() -> None:
    """Objects serialized with the previous module path remain loadable."""
    legacy = b"cMomentRT.phase_function\nIsotropic\n)R."
    assert isinstance(pickle.loads(legacy), Isotropic)
    assert isinstance(pickle.loads(pickle.dumps(Rayleigh())), Rayleigh)
