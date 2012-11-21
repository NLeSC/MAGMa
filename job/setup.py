import os
from setuptools import setup, find_packages
from setuptools.extension import Extension

here = os.path.abspath(os.path.dirname(__file__))
try:
    README = open(os.path.join(here, 'README.rst')).read()
except IOError:
    README = ''


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
    version="1.2",
    license='commercial',
    author='Lars Ridder',
    author_email='lars.ridder@esciencecenter.nl>',
    url='http://www.esciencecenter.nl',
    description='Ms Annotation based on in silico Generated Metabolites',
    long_description=README,
    packages=find_packages(),
    install_requires=['sqlalchemy', 'lxml', 'numpy', 'pp'],
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
)
