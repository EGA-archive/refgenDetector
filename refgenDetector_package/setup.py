from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'refgenDetector'

# Setting up
setup(
    name="refgenDetector",
    version=VERSION,
    author="Mireia Marin Ginestar",
    author_email="<mireia.marin@crg.eu>",
    description=DESCRIPTION,
    install_requires=['argparse', 'pysam', 'gzip'],
    keywords=['python'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix"],
    entry_points={
        'console_scripts': [
            'refgenDetector=refgenDetector.refgenDetector_main:main',
        ],
    },
    packages=find_packages(where='src'),
    package_dir={'': 'src'}

)
