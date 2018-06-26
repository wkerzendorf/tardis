***************************
Updating the reference data
***************************

.. danger::

    UPDATING THE REFERENCE DATA SHOULD ONLY BE DONE IF THE GROUP AGREES THAT THE NEW VERSION
    IS AN IMPROVEMENT.

Generating Plasma Reference
===========================

TARDIS uses py.test to run the tests and generate the reference data. It uses the `plasma` fixture in the
`TestPlasma` class to either run the code or

You can generate Plasma Reference by the following command

.. code-block:: shell

    > pytest -rs tardis/plasma/tests/test_complete_plasmas.py
    --tardis-refdata="/path/to/tardis-refdata/" --generate-reference