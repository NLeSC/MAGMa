from setuptools import setup, find_packages

setup(
    name='Magma',
    version="1.0",
    license='commercial',
    author='Lars Ridder',
    author_email='lars.ridder@esciencecenter.nl>',
    url='http://www.esciencecenter.nl',
    packages=find_packages(),
    install_requires = [ 'sqlalchemy', 'lxml' ],
    dependency_links = [ 'http://www.rdkit.org' ],
    package_data = {
        'magma': ['data/*.smirks', 'script/reactor'],
        },
    entry_points = {
      'console_scripts': [
      	'sygma = magma.script.sygma:main',
        'mscore_mzxml = magma.script.mscore_mzxml:main'
      ],
    }
)
