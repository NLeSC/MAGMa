.. _user:

User management
===============

Users can be added/edited from python prompt.

Start with:

.. code-block:: python

    from transaction import commit
    from sqlalchemy import create_engine
    from magmaweb.user import init_user_db, User
    engine = create_engine('sqlite:///data/users.db')
    init_user_db(engine)


Add user
--------

.. code-block:: python

   user = User(u'stefanv2', u'Stefan Verhoeven',
               u's.verhoeven@esciencecenter.nl', 'mypassword')
   User.add(user)
   commit()

Fetch user
----------

.. code-block:: python

   user = User.by_id('stefanv2')

Update user
-----------

.. code-block:: python

   user.displayname = 'Stefan second account'
   User.add(user)
   commit()

Change password
---------------

.. code-block:: python

   user.password = 'mypw'
   User.add(user)
   commit()

Add job
-------

Jobs started from the web-interface will be registered.
Otherwise the job has to be registered manually.

.. code-block:: python

   from magmaweb.user import JobMeta
   job = JobMeta('00000000-eeee-0000-eeee-000000000000', 'stefanv2',
                 description=u'Same as run[0].description',
                 ms_filename='Same as run[0].ms_filename')
   JobMeta.add(job)
   commit()

Alter owner of job
------------------

.. code-block:: python

   job = JobMeta.by_id('00000000-eeee-0000-eeee-000000000000')
   job.owner = 'someone'
   Job.add(job)
   commit()

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
