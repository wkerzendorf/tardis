"""
Geometry classes for the TARDIS model IO.
"""
from enum import Enum, auto
from dataclasses import dataclass
from typing import Tuple

import numpy as np
from astropy import units as u


class DiscretizationLocation(Enum):
    """
    Enum to define where properties are defined on the mesh.
    """

    CELL_CENTER = auto()
    INTERFACE = auto()


@dataclass
class Radial1DMesh:
    """
    A generic 1D radial mesh defined by spatial coordinates at the interfaces.

    Parameters
    ----------
    radius : astropy.units.Quantity
        The spatial coordinates of the shell interfaces.
    """

    radius: u.Quantity

    @property
    def location(self) -> DiscretizationLocation:
        """
        The location where the primary coordinates (radius) are defined.
        For this mesh, it is always the interfaces.
        """
        return DiscretizationLocation.INTERFACE

    def coordinates(self, location: DiscretizationLocation) -> u.Quantity:
        """
        Get the spatial coordinates for a specific discretization location.

        Parameters
        ----------
        location : DiscretizationLocation
            The requested location (e.g., CELL_CENTER, INTERFACE).

        Returns
        -------
        astropy.units.Quantity
            The coordinates at the requested location.
        """
        if location == DiscretizationLocation.INTERFACE:
            return self.radius
        elif location == DiscretizationLocation.CELL_CENTER:
            return 0.5 * (self.radius[1:] + self.radius[:-1])
        else:
            raise ValueError(f"Location {location} is not supported by Radial1DMesh")


@dataclass
class HomologousRadial1DMesh:
    """
    An intermediate IO mesh defined by velocity interfaces and time.
    Used to bridge legacy/homologous formats to the spatial Radial1DMesh.

    Parameters
    ----------
    velocity : astropy.units.Quantity
        The velocity coordinates of the shell interfaces.
    time_explosion : astropy.units.Quantity
        The time of explosion used to convert velocity to radius.
    """

    velocity: u.Quantity
    time_explosion: u.Quantity

    @property
    def location(self) -> DiscretizationLocation:
        """
        The location where the primary coordinates (velocity) are defined.
        """
        return DiscretizationLocation.INTERFACE

    def to_spatial(self) -> Tuple[Radial1DMesh, u.Quantity]:
        """
        Decompose this homologous mesh into a purely spatial mesh and a velocity field.

        Returns
        -------
        mesh : Radial1DMesh
            The spatial mesh at t = time_explosion.
        velocity : astropy.units.Quantity
            The velocity field defined at the interfaces.
        """
        radius = (self.velocity * self.time_explosion).to(u.cm)
        mesh = Radial1DMesh(radius=radius)
        return mesh, self.velocity
