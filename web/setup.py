import os
import sys

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(here, 'README.rst')).read()
CHANGES = open(os.path.join(here, 'CHANGES.txt')).read()

requires = ['pyramid',
            'SQLAlchemy',
            'pyramid_tm',
            'pyramid_multiauth',
            'pyramid_ipauth',
            'zope.sqlalchemy',
            'colander',
            'waitress',
            'py-bcrypt',
            'Paste'
            ]

tests_require = ['WebTest',
                 'mock'
                 ]

if sys.version_info[:3] < (2, 5, 0):
    requires.append('pysqlite')

setup(name='MAGMaWeb',
      version='0.0',
      description='MAGMaWeb',
      long_description=README + '\n\n' + CHANGES,
      classifiers=["Programming Language :: Python",
                   "Programming Language :: JavaScript",
                   "Intended Audience :: Science/Research",
                   "Environment :: Web Environment",
                   "Natural Language :: English",
                   "Framework :: Pyramid",
                   "Operating System :: OS Independent",
                   "Topic :: Internet :: WWW/HTTP :: WSGI :: Application",
                   "Topic :: Scientific/Engineering :: Chemistry",
                   ],
      author='Stefan Verhoeven',
      author_email='s.verhoeven@esciencecenter.nl',
      url='http://www.emetabolomics.org',
      keywords='web wsgi bfg pylons pyramid cheminformatics',
      packages=find_packages(),
      include_package_data=True,
      zip_safe=False,
      test_suite='magmaweb.tests',
      install_requires=requires,
      tests_require=tests_require,
      entry_points="""\
      [paste.app_factory]
      main = magmaweb:main
      """,
      paster_plugins=['pyramid'],
      )
