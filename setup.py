from setuptools import setup

from cryptic_ecoffs import __version__

setup(
    name='CRyPTIC ECOFF/ECV study ',
    version=__version__,
    author='Philip W Fowler',
    author_email="philip.fowler@ndm.ox.ac.uk",
    description="Data and code to reproduce the figures and tables in the CRyPTIC ECOFF/ECV paper",
    url="https://github.com/oxfordmmm/cryptic-ecoffs",
    packages=['cryptic_ecoffs'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"],
    python_requires='>=3.8',
    install_requires=[
        "numpy",
        "matplotlib",
        "pandas",
        'geopandas',
        'shapely'
    ],
    license="TBD",
    zip_safe=False
)
