from setuptools import setup, find_packages
from nanomgt import version
__version__ = version.__version__
##


setup(
    name='NanoMGT',
    version=__version__,
    packages=find_packages(),
    include_package_data=True,  # Ensure this is set to True
    package_data={
        '': ['*.json'],  # Include all json files from any package
        'nanomgt': ['*.json'],  # Specifically include json files from the nanomgt package if needed
    },
    url='https://github.com/genomicepidemiology/nanomgt',
    license='',
    install_requires=[],
    author='Malte B. Hallgren',
    scripts=['bin/nanomgt'],
    author_email='malhal@food.dtu.dk',
    description='NanoMGT - Nanopore Marker Gene Typing'
)
