from setuptools import setup, find_packages
# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

VERSION = '3.0.0'
DESCRIPTION = 'RefgenDetector'

# Setting up
setup(
    name="RefgenDetector",
    version=VERSION,
    author="Mireia Marin i Ginestar",
    author_email="<mireia.marin@crg.eu>",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=['argparse', 'pysam', 'psutil', 'rich', 'pandas', 'dnspython', 'msgpack'],
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
