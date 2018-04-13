from setuptools import setup, find_packages

setup(
    name = 'primitive_streak_model',
    version = '0.0.1',
    url = 'https://github.com/catohaste/primitive-streak.git',
    author = 'Cato Hastings',
    author_email = 'c.hastings.16@ucl.ac.uk',
    description = 'Description of my package',
    packages = find_packages(),
    install_requires = ['numpy >= 1.11.1', 'matplotlib >= 1.5.1'],
)