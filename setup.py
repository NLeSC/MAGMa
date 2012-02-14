from setuptools import setup, find_packages

setup(
    name='Magma',
    version="1.0",
    license='commercial',
    author='Lars Ridder',
    author_email='lars.ridder@esciencecenter.nl>',
    url='http://www.esciencecenter.nl',
    description='Ms Annotation based on in silico Generated Metabolites',
    packages=find_packages(),
    install_requires=[ 'sqlalchemy', 'zope.sqlalchemy', 'lxml', 'numpy', 'transaction' ],
    dependency_links=[ 'http://www.rdkit.org' ],
    package_data={
        'magma': ['data/*.smirks', 'script/reactor'],
        },
    entry_points={
      'console_scripts': [
        'magma = magma.script:main',
      	'sygma = magma.script.sygma:main',
        'mscore_mzxml = magma.script.mscore_mzxml:main'
      ],
    }
)
