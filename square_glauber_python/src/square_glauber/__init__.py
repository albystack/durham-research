"""Square-grid random-bond dimer model with random-scan heat-bath dynamics."""

from .environment import DimerEnvironment, sample_environment
from .glauber import ChainResult, heatbath_update, run_chain
from .height import face_heights
from .matching import DimerState
from .perfect import CFTPNotCoalesced, MonotoneCFTPSampler, PerfectSample

__all__ = [
    "ChainResult",
    "DimerEnvironment",
    "DimerState",
    "CFTPNotCoalesced",
    "MonotoneCFTPSampler",
    "PerfectSample",
    "face_heights",
    "heatbath_update",
    "run_chain",
    "sample_environment",
]

__version__ = "0.1.0"
