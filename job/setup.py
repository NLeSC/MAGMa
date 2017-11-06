import os
from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy as np

here = os.path.abspath(os.path.dirname(__file__))
try:
    README = open(os.path.join(here, 'README.rst')).read()
except IOError:
    README = ''

# Test if rdkit is present with INCHI support
try:
    from rdkit.Chem.inchi import INCHI_AVAILABLE
    if not INCHI_AVAILABLE:
        raise Exception('RDKit with INCHI support is required')
except ImportError:
    raise Exception('RDKit with INCHI support is required')

# Test if pp is present
try:
    import pp
except ImportError:
    raise Exception('Parallel Python (pp) is required')

# Only use Cython if it is available, else just use the pre-generated files
try:
    from Cython.Distutils import build_ext
    source_ext = '.pyx'
    cmdclass = {'build_ext': build_ext}
except ImportError:
    # If missing can be created with 'cython magma/fragmentation_cy.pyx'
    source_ext = '.c'
    cmdclass = {}

ext_modules = [Extension('magma.fragmentation_cy',
                         ['magma/fragmentation_cy' + source_ext])]

setup(
    name='Magma',
    version='1.3',
    license='commercial',
    author='Lars Ridder',
    author_email='lars.ridder@esciencecenter.nl>',
    url='http://www.esciencecenter.nl',
    description='Ms Annotation based on in silico Generated Metabolites',
    long_description=README,
    classifiers=["Intended Audience :: Science/Research",
                   "Environment :: Console",
                   "Natural Language :: English",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Chemistry",
                   ],
    packages=find_packages(),
    install_requires=['sqlalchemy', 'lxml', 'numpy', 'requests', 'macauthlib', 'mock', 'nose', 'coverage'],
    package_data={
        'magma': ['data/*.smirks', 'script/reactor'],
        },
    entry_points={
      'console_scripts': [
        'magma = magma.script:main',
      ],
    },
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    include_dirs=[np.get_include()],
)
