"""Magma pyramid web app"""
from pyramid.config import Configurator


def main(global_config, **settings):
    """This function returns the Magma WSGI application."""
    # pylint: disable=W0613
    config = Configurator(settings=settings)
    config.include('magmaweb.config.configure')
    return config.make_wsgi_app()
