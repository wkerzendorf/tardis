.. _virtual_packets:

****************************************
Spectrum Generation with Virtual Packets
****************************************

The main purpose of TARDIS is the generation of synthetic spectra. 

Spectrum Generation Methods
============================

TARDIS provides three methods for generating synthetic spectra:

1. **Real Packet Spectrum**: Records the properties of all escaping Monte Carlo
   packets and bins their contributions in frequency space. This direct approach
   naturally suffers from Monte Carlo noise.

2. **Formal Integral Spectrum**: Uses the converged radiation field and plasma
   state to analytically compute the emergent spectrum via the formal solution
   of the radiative transfer equation. See :ref:`formal_integral` for details.

3. **Virtual Packet Spectrum**: Employs a variance reduction technique where
   additional "virtual" packets are traced from real packet spawn events. This
   method provides smoother spectra with better signal-to-noise while maintaining
   physical accuracy.

This document describes the virtual packet method in detail.

Monte Carlo Noise and Variance Reduction
=========================================

One implementation follows the obvious approach of recording the properties
of all escaping Monte Carlo packets and binning their contributions in
frequency (or wavelength) space. This "real packet" spectrum will naturally
suffer from Monte Carlo noise, and if one tries to improve its signal-to-noise
ratio, one immediately encounters a fundamental characteristic of Monte Carlo
approaches. Since Monte Carlo processes are typically governed by Poisson
statistics, the level of stochastic fluctuations decreases only as :math:`\propto
N^{-\frac{1}{2}}`, with :math:`N` denoting the number of Monte Carlo
experiments. Thus, to decrease the noise by an order of magnitude, 100 times
more experiments have to be performed. In the case of the real packet spectrum,
this translates into using 100 times more packets, which would increase the
runtime by about the same factor.

.. note::

    More details about Monte Carlo errors and noise behaviour may be found in
    the standard literature, for example in :cite:`Kalos2008`.

It is difficult to avoid this fundamental behaviour of Monte Carlo techniques.
However, sophisticated Monte Carlo techniques exist, which make better use of
the computational resources. One such approach, which achieves a better noise
behaviour for the same computational costs, is implemented in TARDIS. It relies
on the concept of so-called "virtual packets" and goes back to the works by
:cite:`Long2002` and :cite:`Sim2010`.

Virtual Packets
===============

The virtual packet scheme is best explained by first detailing how it works and
then illustrating the physical reasoning for introducing this scheme. More
information about this scheme may be found in :cite:`Kerzendorf2014`.

Virtual Packet Procedure
------------------------

In the virtual packet scheme, a new population of Monte Carlo packets is
introduced. Every time a real packet is launched or performs a physical
interaction, a pre-defined number of virtual packets, :math:`N_v`, are
generated. These propagation of these "virtual packets" is followed in a
similar fashion to the real ones with the important distinction that their
trajectory is never changed. However, the optical depth the virtual packet
accumulates during its propagation to the ejecta surface due to electron
scattering and atomic line interactions is tracked. Once the virtual packet
escapes, its contribution to the spectrum is then weighted by this total
optical depth. In particular, it contributes with the

.. math::

    \Delta L_{\nu} = \varepsilon_{\nu} \exp(-\tau) \frac{1}{\Delta t \Delta \nu}

to the emergent luminosity in the frequency interval :math:`[\nu, \nu + \Delta
\nu]`. Here, :math:`\Delta t` denotes the physical duration of the simulation
step (the same duration which is used during the initialization process at the
photosphere, see :ref:`Propagation <propagation>`), and :math:`\varepsilon` is
the energy of the virtual packet when it was generated.

.. note::

    TARDIS is a time-independent radiative transfer scheme. Thus, :math:`\Delta
    t` should be interpreted as the physical duration of the emission process
    at the photosphere.


The initialization process for virtual packets is slightly different from real
ones. For example, whenever a real packet is emitted by the photosphere,
:math:`N_v` virtual packets are spawned, as well. The propagation direction of
these virtual packets is assigned uniformly. However, since :math:`N_v` is
typically small, an unequal sampling of the solid angle is avoided by selecting
:math:`N_v` equal :math:`\mu` bins and sampling the direction uniformly within
these bins. Since the emitted radiation field has a different angular
dependence (represented by the non-uniform sampling rule for real packets,
:math:`\mu = \sqrt{z}`), the energy of each virtual packet is weighted accordingly

.. math::

    \varepsilon_v = \varepsilon \frac{2 \mu}{N_v}.

Here, :math:`\varepsilon` is the energy of the real packet that spawned the
virtual ones. If virtual packets are generated during a real packet interaction
at the location :math:`r`, their propagation direction is sampled uniformly
from the interval :math:`[\mu_{\mathrm{min}}, 1]`. 

.. math::

    \mu_{\mathrm{min}} = - \sqrt{1 - \left(\frac{R_{\mathrm{phot}}}{r}\right)^2}

Setting this lower limit avoids virtual packet trajectories intercepting the
photosphere at :math:`R_{\mathrm{phot}}` (in which case the virtual packet
could not contribute to the emergent spectrum and computational resources would
have been wasted on this packet). The amount of radiation being backscattered
towards the photosphere is accounted for by modifying the energy of the virtual
packets

.. math::

    \varepsilon_v = \varepsilon \frac{1 - \mu_{\mathrm{min}}}{2 N_v}.


Interpretation
--------------

The basic idea underlying the virtual packet scheme may be illustrated by
considering the formal solution of the time-independent radiative transfer
problem :

.. math::

    I(R, \mu, \nu) = I(R_{\mathrm{phot}}, \mu, \nu) \exp(-\tau(s_0)) +
    \int_0^{s_0} \eta(R - \mu s, \mu, \nu) \exp(-\tau(s)) \mathrm{d}s

This formulation of the formal solution is valid for the supernova ejecta problem and
involves the location of the photosphere, the radius of the ejecta surface
:math:`R` and the packet trajectory :math:`s`. Here, the optical depth
:math:`\tau(s)` measures the optical depth from :math:`s` to the ejecta surface.
For more details see :ref:`montecarlo_basics`.

Essentially, the virtual packets solve this formal solution equation along a
large number of directional rays. In particular, the virtual packets spawned at
the photosphere solve the first part of the formal solution, namely by
determining which fraction of the photospheric radiation field remains at the
surface of the ejecta. The virtual packets which are generated whenever a real
packet interacts, account for the second part of the formal solution. In this
interpretation, the purpose of the real packet population is simply to "sample"
the emissivity of the medium.

This outline of the virtual packet scheme is concluded with a remark about its
benefits. The advantages of using a combination of real and virtual packets
compared to calculation based purely on real packets lies in lower
computational costs which are associated with solving the propagation of
virtual packets. These always propagate along a straight line, whereas real
packets may be deflected multiple times, thus making the determination of the
entire propagation path more expensive.


Virtual Packet Generation Workflow
===================================

Implementation via Tracker-Based Postprocessing
------------------------------------------------

In TARDIS, virtual packets are generated via **postprocessing** after the Monte
Carlo transport simulation completes. This approach decouples spectrum generation
from the transport calculation and enables performance optimizations.

The workflow consists of two phases:

**Phase 1: Monte Carlo Transport with Tracking**

When virtual packets are requested (``no_of_virtual_packets > 0``), TARDIS
automatically enables full real packet tracking. During transport, every real
packet interaction that could spawn virtual packets is logged to the tracker::

    tracker_full_df

This tracker records "spawn events" including:

- Initial packet emission from the photosphere
- Line interactions (resonant scattering events)
- Electron scattering (ESCATTER) events

Each spawn event contains the packet state (position, direction, energy, frequency)
at the moment virtual packets would be generated.

**Phase 2: Virtual Packet Generation from Tracker**

After the simulation completes, virtual packets are generated on-demand by
calling::

    sim.generate_virtual_spectrum()

This method:

1. Extracts spawn events from ``tracker_full_df``
2. For each spawn event, generates :math:`N_v` virtual packets
3. Traces each virtual packet through the ejecta (accumulating optical depth)
4. Bins the emergent virtual packet energies to create the spectrum

This postprocessing approach uses numba parallelization to achieve 10-30%
performance improvements compared to inline generation.

Usage Example
-------------

.. code-block:: python

    from tardis import run_tardis
    from tardis.io.configuration.config_reader import Configuration

    # Load configuration with virtual packets enabled
    config = Configuration.from_yaml("tardis_example.yml")
    # config.montecarlo.no_of_virtual_packets = 10  # if not in YAML

    # Phase 1: Run Monte Carlo transport (tracking enabled automatically)
    sim = run_tardis(config)

    # Phase 2: Generate virtual packets from tracker
    sim.generate_virtual_spectrum()

    # Access the three spectrum types
    spectrum_real = sim.spectrum_solver.spectrum_real_packets
    spectrum_integrated = sim.spectrum_solver.spectrum_integrated  
    spectrum_virtual = sim.spectrum_solver.spectrum_virtual_packets

Benefits of Postprocessing
---------------------------

The tracker-based postprocessing approach provides several advantages:

- **Performance**: Numba-parallelized generation is faster than inline tracking
- **Flexibility**: Generate multiple spectrum variations from a single simulation
- **Memory efficiency**: Tracker is more compact than storing all virtual packet data
- **Decoupling**: Spectrum generation can be optimized independently from transport

.. note::

    The physical behavior of virtual packets (angular sampling, energy weighting,
    optical depth accumulation) remains identical to previous versions. Only the
    *timing* of generation has changed from inline to postprocessing.

For implementation details, see :py:class:`~tardis.spectrum.virtual_packet_solver.VirtualPacketSolver`
and :py:class:`~tardis.spectrum.virtual_packet_state.VirtualPacketState`.

