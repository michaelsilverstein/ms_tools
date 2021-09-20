from setuptools import setup
from ms_tools import __version__

dependencies = ['numpy', 'pandas']

with open('README.md') as fh:
    long_description = fh.read()

setup(
    name='multiloan',
    version=__version__,
    packages=['ms_tools'],
    install_requires=dependencies,
    url='',
    license='MIT',
    author='michaelsilverstein',
    author_email='michael.silverstein4@gmail.com',
    description='Functions that I commonly use',
    long_description = long_description,
    long_description_content_type = "text/markdown",
)