.. _faq:

**************************
Frequently Asked Questions
**************************

**This reference section provides quick answers to common questions about TARDIS.**

This FAQ serves as a reference for finding solutions to frequently encountered issues. For step-by-step learning, see our :doc:`tutorials`. For specific problem-solving guides, check our :doc:`how_to_guides`.

Overview
--------

Here you'll find answers to common questions organized by topic.

- :ref:`faq-usage`

     - :ref:`faq-usage-memory`

- :ref:`faq-spectrum`

     - :ref:`faq-spectrum-vpacket-none`
     - :ref:`faq-spectrum-which-method`
     - :ref:`faq-spectrum-visualization-fails`


.. _faq-usage:

Usage
-----

.. _faq-usage-memory:

My simulation seems to consume excessive amounts of memory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


High memory usage can have many reasons. Both your model
and certain simulation settings can increase the memory
usage significantly

1. Enabling ``track_rpacket: true`` will take up substantial
   amounts of memory, in particular for higher packet counts.
   Consider turning this feature off or only use it with a
   small amount of :term:`packets`.
2. Both the number of shells and the number of :term:`packets`
   increase the memory requirements. Consider using only
   as many shells/ :term:`packets` as are required to converge
   on a result.

.. _faq-spectrum:

Spectrum Generation
-------------------

.. _faq-spectrum-vpacket-none:

Virtual packet spectrum is None or AttributeError occurs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you're trying to access virtual packet spectrum or use visualization tools and encounter:

.. code-block:: python

    spectrum = sim.spectrum_solver.spectrum_virtual_packets
    # AttributeError: 'NoneType' object has no attribute 'plot'

**Solution:** Virtual packets now require an explicit postprocessing step after simulation:

.. code-block:: python

    sim = run_tardis("config.yml")
    sim.generate_virtual_spectrum()  # Add this line
    spectrum = sim.spectrum_solver.spectrum_virtual_packets

This applies to workflows as well:

.. code-block:: python

    workflow = StandardWorkflow()
    workflow.run(config="config.yml")
    workflow.generate_virtual_spectrum()  # Add this line

See the `Virtual Packet Migration Guide <conceptual_changelog/2026_virtual_packet_refactor.md>`_ for complete migration instructions.

.. _faq-spectrum-which-method:

Which spectrum generation method should I use?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TARDIS offers three spectrum generation methods:

1. **Real Packet Spectrum** (``spectrum_real_packets``): Direct Monte Carlo output
   
   - **When to use**: Quick checks, debugging
   - **Pros**: No additional computation
   - **Cons**: Noisy, requires many packets for smooth spectrum

2. **Formal Integral Spectrum** (``spectrum_integrated``): Analytical post-processing
   
   - **When to use**: High-quality spectra, publication plots
   - **Pros**: Smooth, physically consistent with converged radiation field
   - **Cons**: Computationally expensive

3. **Virtual Packet Spectrum** (``spectrum_virtual_packets``): Variance reduction technique
   
   - **When to use**: Balance between smoothness and computational cost
   - **Pros**: Moderate computational cost, smoother than real packets, enables SDEC/LIV analysis
   - **Cons**: Requires postprocessing call

**Recommendation:** Use virtual packets (method 3) for most analysis work, and formal integral for final publication-quality spectra.

.. _faq-spectrum-visualization-fails:

Visualization tools (SDEC, LIV plots) fail to work
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If SDECPlotter, LIVPlotter, or widget generation fails with errors like:

.. code-block:: text

    Error: virtual_packet_state is None

**Solution:** Generate virtual packets before creating visualizations:

.. code-block:: python

    from tardis.visualization import SDECPlotter
    from tardis import run_tardis
    
    sim = run_tardis("config.yml")
    sim.generate_virtual_spectrum()  # Required for visualization tools
    
    plotter = SDECPlotter.from_simulation(sim)
    plotter.generate_plot_mpl()

All TARDIS visualization tools (SDEC, LIV, widgets) require virtual packets and thus need the ``generate_virtual_spectrum()`` call.

See the :doc:`../analyzing_tardis/visualization/index` for visualization tutorials.
