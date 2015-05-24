Helium Non-Thermal Treatment Approximation (ab-NLTE)
------------------------------------------

This module will be introduced as part of the new plasma structure. It uses a simple approximation to account for changes in the state of helium due to non-thermal collisions with fast electrons. The detailed physics behind this will be published in an upcoming astronomy paper.

Unlike the other elements in Tardis, the state of helium is primarily determined by the rate of collisions with fast electrons generated as a result of radioactive decay, rather than by the rate of excitation or ionisation by photons. This is because helium takes an unusually large amount of energy to excite or ionise. Generally the rate of fast electron collisions and the resulting state of helium would be calculated using numerical simulations, but we have found a much simpler approximation.

The approximation makes a few assumptions that are all very well supported by numerical simulations. The approximation first assumes that the vast majority of helium ions exist in the He II ground state, and that the He I ground state has a population of zero. The He II excited states and He III ground state are calculated relative to the He II ground state population using the typical equations of dilute LTE. The He I excited states are then calculated relative to the He II ground state population as well, as though they were in dilute equilibrium with the He II ground state rather than their own ground state. 

.. math::

The He I excited level populations are calculated relative to the He II ground state population as follows: 

\frac{n_{i,j,k}}{n_{i,j+1,0}n_{e}} = \frac{1}{W}\frac{g_{i},j}{2g_{i,j+1,0}}\left(\frac{h^{2}}{2\pi m_{e}kT_{r}}\right)^{\frac{3}{2}}e^{\frac{\chi_{i,j}-\epsilon_{i,j,k}}{kT_{r}}}

n: level population
i: atomic number
j: ion number
k: level number
n_{e}: electron density
W: dilution factor
g: statistical weight
T_{r}: radiation temperature
\chi: ionisation energy
\epsilon: excitation energy

The He II excited state and He III ground state populations are calculated using the same equations as are behind the dilute LTE excitation and nebular ionisation plasma modes. However, values are calculated relative to the He II ground state population rather than to the ion population as a whole, using the fact that:

\frac{N_{i,j}}{Z_{i,j}} = \frac{n_{i,j,0}}{g_{i,j,0}}
