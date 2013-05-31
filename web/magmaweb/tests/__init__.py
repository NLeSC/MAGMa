"""Tests for MAGMa web"""
from unittest import TestCase
from mock import patch
from magmaweb import main


class TestMain(TestCase):

    def testIt(self):
        settings = {'restricted': True}

        with patch('magmaweb.Configurator', spec=True) as Configurator:

            app = main({}, **settings)

            Configurator.assert_called_with(settings=settings)
            config = Configurator.return_value
            config.include.assert_called_with('magmaweb.config.configure')
            config.make_wsgi_app.assert_called_with()