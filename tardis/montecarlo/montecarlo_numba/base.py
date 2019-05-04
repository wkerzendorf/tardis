from numba import prange, njit, config
import numpy as np
from tardis.montecarlo.montecarlo_numba.rpacket import RPacket, PacketStatus
from tardis.montecarlo.montecarlo_numba.storage_model import StorageModel, initialize_storage_model
from tardis.montecarlo.montecarlo_numba.single_packet_loop import single_packet_loop
from tardis.montecarlo.montecarlo_numba import njit_dict

#config.THREADING_LAYER = 'threadsafe'
#config.DEBUG_ARRAY_OPT=1

def montecarlo_radial1d(model, plasma, runner):
    storage_model = initialize_storage_model(model, plasma, runner)
    montecarlo_main_loop(storage_model)

@njit(**njit_dict, nogil=True)
def montecarlo_main_loop(storage_model):
    """
    This is the main loop of the MonteCarlo routine that generates packets 
    and sends them through the ejecta. 

    Parameters
    ----------
    storage_model : [type]
        [description]
    """
    output_nus = np.empty_like(storage_model.output_nus)
    output_energies = np.empty_like(storage_model.output_energies)
    for i in prange(storage_model.no_of_packets):
        r_packet = RPacket(storage_model.r_inner[0],
                           storage_model.packet_mus[i],
                           storage_model.packet_nus[i],
                           storage_model.packet_energies[i])
        single_packet_loop(storage_model, r_packet)
        output_nus[i] = r_packet.nu

        if r_packet.status == PacketStatus.REABSORBED:
            output_energies[i] = -r_packet.energy
        elif r_packet.status == PacketStatus.EMITTED:
            output_energies[i] = r_packet.energy
    storage_model.output_energies[:] = output_energies[:]
    storage_model.output_nus[:] = output_nus[:]

