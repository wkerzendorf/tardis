Integration Tests
*****************

The integration tests are more advanced tests that take several hours to days to run. There is a custom setup for TARDIS
to work with them.

Running the integration tests
=============================

These tests require reference files against which the results of the various
tardis runs are tested. So you first need to either download the current
reference files (`here <https://github.com/tardis-sn/tardis-refdata>`_)
or generate new ones.

Both of of these require a configuration file for the integration tests:

.. literalinclude:: yml_files/integration.yml
    :language: yaml

Inside the atomic data directory there needs to be atomic data for each of
the setups that are provided in the ``test_integration`` folder.
If no references are given the first step is to generate them.
The ``--less-packets`` option is useful for debugging purposes and will just
use very few packets to generate the references and thus make the process much
faster - THIS IS ONLY FOR DEBUGGING PURPOSES. The ``-s`` option ensures that
TARDIS prints out the progress:

.. code-block::

    > python setup.py test --args="--integration=integration.yml -m integration
    --generate-reference --less-packets"

To run the test after having run the ``--generate-references`` all that is
needed is:

.. code-block::

    > python setup.py test --args="--integration=integration.yml -m integration
    --less-packets" --remote-data


Setting up the Dokuwiki report
==============================

A normal dokuwiki installation is performed on the required server. Before the
connection works one is requires to set the option remote access in the
settings. If this is not done the ``dokuwiki`` python plugin will not connect
with the warning ``DokuWikiError: syntax error: line 1, column 0``. One also
has to enable this for users (``remoteuser`` option) otherwise the error:
``ProtocolError for xmlrpc.php?p=xxxxxx&u=tardistester: 403 Forbidden``
will appear.

Another important configuration option is to enable embedded html ``htmlok``
otherwise it won't show nice html page reports.

Finally, one has to call the `python setup.py test` with the ``--remote-data``
option to allow posting to an external dokuwiki server.
