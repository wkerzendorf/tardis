0.9.2 (unreleased)
------------------

Bugfixes
^^^^^^^^

- "print" replaced with logger for classical nebular
- logger statement added for coronal approximation
- warning added to documentation since plasma is out of date (temp
  solution only) #108
- bug fix in plasma_array, wrong indentation on line 587
  (calculate_nlte_level_populations)
- simplification of nlte rates - treatment of stimulated emission no
  longer handled as negative absorption (correction factor removed
  from calculate_nlte_level_populations / reformulated terms)


0.9.1 (2014-02-03)
------------------

New Features
^^^^^^^^^^^^

- bugfix release of TARDIS
- updated example section in the documentation
- now with working test coverage setup (https://coveralls.io/r/tardis-sn/tardis)


Bugfixes
^^^^^^^^

- missing ez_setup.py and setuptools_bootstrap.py included now
- reading of ascii density and abundance profiles
- several fixes in the documentation


