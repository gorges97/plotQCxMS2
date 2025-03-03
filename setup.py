# -*- coding: utf-8 -*-


from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='plotQCxMS2',
    version='0.1.0',
    description='Plotting tool for QCxMS2',
    long_description=readme,
    author='Johannes Gorges',
    author_email='gorges@thch.uni-bonn',
    url='https://github.com/gorges97/plotQCxMS2',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

