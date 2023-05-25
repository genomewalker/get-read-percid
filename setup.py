from setuptools import setup
import versioneer

requirements = [
    "pandas>=1.4.2",
    "tqdm>=4.64.1",
    "pysam>=0.17.0",
    "numpy>=1.21.2",
    "biopython>=1.79",
]

setup(
    setup_requires=[
        # Setuptools 18.0 properly handles Cython extensions.
        "setuptools>=39.1.0",
        "Cython>=0.29.24",
    ],
    name="get-read-percid",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A simple tool to extract references from BAM files and get read statistics",
    license="MIT",
    author="Antonio Fernandez-Guerra",
    author_email="antonio@metagenomics.eu",
    url="https://github.com/genomewalker/get-read-percid",
    packages=["get_read_percid"],
    entry_points={"console_scripts": ["getRPercId=get_read_percid.__main__:main"]},
    install_requires=requirements,
    keywords="get-read-percid",
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
