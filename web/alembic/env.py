from __future__ import with_statement
from alembic import context
from sqlalchemy import create_engine, pool
from logging.config import fileConfig

# this is the Alembic Config object, which provides
# access to the values within the .ini file in use.
config = context.config

# Interpret the config file for Python logging.
# This line sets up loggers basically.
fileConfig(config.config_file_name)

# add your model's MetaData object here
# for 'autogenerate' support
# from myapp import mymodel
# target_metadata = mymodel.Base.metadata
#target_metadata = None
from magmaweb.models import Base
target_metadata = Base.metadata

# other values from the config, defined by the needs of env.py,
# can be acquired:
# my_important_option = config.get_main_option("my_important_option")
# ... etc.

def run_migrations_offline():
    """Run migrations in 'offline' mode.

    This configures the context with just a URL
    and not an Engine, though an Engine is acceptable
    here as well.  By skipping the Engine creation
    we don't even need a DBAPI to be available.

    Calls to context.execute() here emit the given string to the
    script output.

    """
    jobid = context.get_x_argument(as_dictionary=True).get("jobid")
    url = config.get_main_option("jobdb").format(jobid=jobid)
    context.configure(url=url)

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online():
    """Run migrations in 'online' mode.

    In this scenario we need to create an Engine
    and associate a connection with the context.

    """
    jobid = context.get_x_argument(as_dictionary=True).get("jobid")
    jobdb = config.get_main_option("jobdb").format()
    url = jobdb.format(jobid=jobid)
    engine = create_engine(
                url,
                poolclass=pool.NullPool)

    connection = engine.connect()
    context.configure(
                connection=connection,
                target_metadata=target_metadata
                )

    try:
        with context.begin_transaction():
            context.run_migrations()
    finally:
        connection.close()

if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()

