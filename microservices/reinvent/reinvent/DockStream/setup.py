from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))


setup(
    name="DockStream",
    version="1.0.0",
    description=(
        "A comprehensive Python wrapper for a selection of molecule docking backends"
    ),
    license="Apache License 2.0",
    url="https://github.com/MolecularAI/DockStream",
    python_requires=">=3.6.1",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.6",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Libraries :: Application Frameworks"
    ],
    keywords="docking, molecular embedding, bio-physics",
    packages=find_packages(exclude=['docs', 'tests*']),
    include_package_data=True,
    author="Christian Margreitter",
    author_email="christian.margreitter@astrazeneca.com"
)
