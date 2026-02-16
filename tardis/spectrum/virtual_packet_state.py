"""Virtual packet state for postprocessing-generated virtual packets."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from astropy import units as u

if TYPE_CHECKING:
    import numpy as np

    from tardis.transport.montecarlo.packets.packet_collections import (
        VPacketCollection,
    )


@dataclass
class VirtualPacketState:
    """
    Holds consolidated virtual packet data generated in postprocessing.

    This class follows the State/Solver pattern used in TARDIS. It stores
    the results of postprocessing virtual packet generation from tracker data.

    Attributes
    ----------
    vpacket_collection : VPacketCollection
        Consolidated collection of all virtual packets with last_interaction metadata.
    luminosity : np.ndarray
        Energy histogram converted to luminosity [erg/s], binned to spectrum frequency grid.
    count : int
        Total number of virtual packets generated.

    Notes
    -----
    This class provides properties for backward compatibility with visualization
    tools (SDEC, LIV plots) that previously accessed vpacket data from
    MonteCarloTransportState.
    """

    vpacket_collection: VPacketCollection
    luminosity: np.ndarray
    count: int

    @property
    def nus(self) -> u.Quantity:
        """
        Virtual packet frequencies.

        Returns
        -------
        astropy.units.Quantity
            Frequencies [Hz].
        """
        return u.Quantity(self.vpacket_collection.nus[: self.count], u.Hz)

    @property
    def energies(self) -> u.Quantity:
        """
        Virtual packet energies.

        Returns
        -------
        astropy.units.Quantity
            Energies [erg].
        """
        return u.Quantity(self.vpacket_collection.energies[: self.count], u.erg)

    @property
    def initial_rs(self) -> u.Quantity:
        """
        Virtual packet initial radii.

        Returns
        -------
        astropy.units.Quantity
            Radii [cm].
        """
        return u.Quantity(
            self.vpacket_collection.initial_rs[: self.count], u.cm
        )

    @property
    def initial_mus(self) -> u.Quantity:
        """
        Virtual packet initial direction cosines.

        Returns
        -------
        astropy.units.Quantity
            Direction cosines (dimensionless).
        """
        return u.Quantity(
            self.vpacket_collection.initial_mus[: self.count],
            u.dimensionless_unscaled,
        )

    @property
    def last_interaction_type(self) -> np.ndarray:
        """
        Last interaction type for each virtual packet's parent rpacket.

        Returns
        -------
        np.ndarray
            Interaction type codes.
        """
        return self.vpacket_collection.last_interaction_type[: self.count]

    @property
    def last_interaction_in_nu(self) -> u.Quantity:
        """
        Frequency before last interaction of parent rpacket.

        Returns
        -------
        astropy.units.Quantity
            Frequencies [Hz].
        """
        return u.Quantity(
            self.vpacket_collection.last_interaction_in_nu[: self.count], u.Hz
        )

    @property
    def last_interaction_in_r(self) -> u.Quantity:
        """
        Radius at last interaction of parent rpacket.

        Returns
        -------
        astropy.units.Quantity
            Radii [cm].
        """
        return u.Quantity(
            self.vpacket_collection.last_interaction_in_r[: self.count], u.cm
        )

    @property
    def last_interaction_in_id(self) -> np.ndarray:
        """
        Inbound atomic level ID for last line interaction of parent rpacket.

        Returns
        -------
        np.ndarray
            Level IDs (-1 for non-line interactions).
        """
        return self.vpacket_collection.last_interaction_in_id[: self.count]

    @property
    def last_interaction_out_id(self) -> np.ndarray:
        """
        Outbound atomic level ID for last line interaction of parent rpacket.

        Returns
        -------
        np.ndarray
            Level IDs (-1 for non-line interactions).
        """
        return self.vpacket_collection.last_interaction_out_id[: self.count]

    @property
    def last_interaction_shell_id(self) -> np.ndarray:
        """
        Shell ID where last interaction of parent rpacket occurred.

        Returns
        -------
        np.ndarray
            Shell IDs.
        """
        return self.vpacket_collection.last_interaction_shell_id[: self.count]
