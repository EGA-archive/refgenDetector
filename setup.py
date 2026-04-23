from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

VERSION = '3.0.1'
DESCRIPTION = 'RefgenDetector'

setup(
    name="RefgenDetector",
    version=VERSION,
    author="Mireia Marin i Ginestar",
    author_email="<mireia.marin@crg.eu>",
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=['pysam', 'psutil', 'rich', 'pandas', 'dnspython', 'msgpack', 'numpy'],
    keywords=['python'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
    ],
    entry_points={
        'console_scripts': [
            'refgenDetector=refgenDetector.refgenDetector_main:main',
        ],
    },
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    package_data={"refgenDetector": ["*.msgpack"]}
)
