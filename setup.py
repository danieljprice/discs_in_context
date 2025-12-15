"""
Setup script for discs_in_context package.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="discs-in-context",
    version="0.1.0",
    author="",
    description="A package for plotting extinction maps of star-forming regions with protostars and protoplanetary discs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'discs_in_context': ['data/*.csv'],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
)

