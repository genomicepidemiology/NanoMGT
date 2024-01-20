from setuptools import setup, find_packages

setup(
    name='NanoMGT',
    version='1.0.0',
    packages=find_packages(),
    data_files=[],
    include_package_data=True,
    url='https://github.com/genomicepidemiology/nanomgt',
    license='',
    install_requires=(),
    author='Malte B. Hallgren',
    scripts=['bin/nanomgt'],
    author_email='malhal@food.dtu.dk',
    description='NanoMGT - Nanopore Marker Gene Typing'
)