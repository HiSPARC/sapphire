from setuptools import setup, find_packages
setup(
    name = "sapphire",
    version = "0.9.2b",
    packages = find_packages(),
    url = "http://github.com/hisparc/sapphire/",
    author = "David Fokkema",
    author_email = "davidf@nikhef.nl",
    description = "A framework for the HiSPARC experiment",

    install_requires = ['numpy', 'scipy', 'tables', 'matplotlib',
                        'progressbar', 'mock'],
)
