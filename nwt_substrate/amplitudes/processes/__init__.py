"""Pre-built substrate-algebra amplitudes for canonical SM processes."""

from . import compton
from . import eemumu
from . import moller
from . import bhabha
from . import muon_decay
from . import qcd_qqbar
from . import qcd_qq
from . import qcd_gg

__all__ = ["compton", "eemumu", "moller", "bhabha", "muon_decay",
           "qcd_qqbar", "qcd_qq", "qcd_gg"]
