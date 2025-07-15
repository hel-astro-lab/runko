from pyrunko._runko_next import *
from .Configuration import Configuration
from .TileGrid import TileGrid
import pyrunko.emf2.threeD as emf
from .Simulation import Simulation
from pyrunko.tools import comm_mode


def on_main_rank() -> bool:
    """
    Checks if the caller is "main" rank.

    There is no other guarantees other that there is only one main rank.
    """

    try:
        from mpi4py import MPI
        return MPI.COMM_WORLD.Get_rank() == 0
    except:
        # If importing mpi4py fails, assume that there is no MPI
        # which means that there has to be only one rank.
        return True

