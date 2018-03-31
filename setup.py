from setuptools import find_packages, setup

# set version number and write to sapphire/version.py
version = '1.6.0'

version_py = """\
# Created by setup.py. Do not edit.
__version__ = "{version}"
"""
with open('sapphire/version.py', 'w') as f:
    f.write(version_py.format(version=version))


setup(name='hisparc-sapphire',
      version=version,
      packages=find_packages(),
      url='http://github.com/hisparc/sapphire/',
      bugtrack_url='http://github.com/HiSPARC/sapphire/issues',
      license='GPLv3',
      author='David Fokkema, Arne de Laat, and others',
      author_email='davidf@nikhef.nl, adelaat@nikhef.nl',
      maintainer='HiSPARC',
      maintainer_email='beheer@hisparc.nl',
      description='A framework for the HiSPARC experiment',
      long_description=open('README.rst').read(),
      keywords=['HiSPARC', 'Nikhef', 'cosmic rays'],
      classifiers=[
          'Intended Audience :: Science/Research',
          'Intended Audience :: Education',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering',
          'Topic :: Education',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'],
      scripts=['sapphire/corsika/generate_corsika_overview',
               'sapphire/corsika/qsub_corsika',
               'sapphire/corsika/qsub_store_corsika_data',
               'sapphire/corsika/store_corsika_data',
               'sapphire/data/update_local_data',
               'sapphire/data/extend_local_data'],
      package_data={'sapphire': ['data/*.json',
                                 'data/*/*.json',
                                 'data/current/*.tsv',
                                 'data/detector_timing_offsets/*.tsv',
                                 'data/electronics/*.tsv',
                                 'data/gps/*.tsv',
                                 'data/layout/*.tsv',
                                 'data/station_timing_offsets/*/*.tsv',
                                 'data/trigger/*.tsv',
                                 'data/voltage/*.tsv',
                                 'corsika/LICENSE',
                                 'tests/test_data/*.h5',
                                 'tests/test_data/*.tsv',
                                 'tests/test_data/*.dat',
                                 'tests/analysis/test_data/*.h5',
                                 'tests/corsika/test_data/*.h5',
                                 'tests/corsika/test_data/*/DAT000000',
                                 'tests/corsika/test_data/*/*.h5',
                                 'tests/simulations/test_data/*.h5']},
      install_requires=['numpy', 'scipy', 'tables>=3.3.0',
                        'progressbar2>=3.7.0', 'lazy', 'mock', 'six'],
      extras_require={
          'dev': ['Sphinx', 'flake8', 'pep8-naming', 'coverage', 'flake8-isort'],
          'astropy': ["astropy"]},
      test_suite="sapphire.tests",)
