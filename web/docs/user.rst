.. _user:

User & job management
=====================

Users and jobs can be added/edited/deleted using the `magma-web` script.

Anonymous and restricted mode
=============================

To make the MAGMa application usable by everyone,
MAGMa can be configured to run in anonymous and restricted mode.

Also for certain journals it is required that you can use the application without creating an account.

Anonymous mode
--------------

To enable put in the *.ini config file:

.. code-block:: ini

  auto_register = true

If you go to the application you will not be asked to login.
You still have a workspace and make results public if you want.

Restricted mode
---------------

To enable put in the *.ini config file:

.. code-block:: ini

  restricted = true

This applies restrictions to the calculations to make sure calculations run fast.
