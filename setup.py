import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hlaquant",
    version="0.0.1",
    author="Austin Crinklaw",
    author_email="acrinklaw@lji.org",
    description="HLAQuant - Get HLA allele specific expression",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/acrinklaw/hlaquant",
    packages=setuptools.find_packages(),
    package_data={'HLAQuant': ['data/blastdb/*', 'data/hla_nuc.fasta']},
    install_requires=[
        'pandas',
        'biopython',],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: Unix",
    ],
)