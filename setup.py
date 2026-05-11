from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()
VERSION = '3.0.5'
DESCRIPTION = 'RefgenDetector'


# ── Post-install hook ─────────────────────────────────────────────────────────

def run_post_install():
    try:
        from refgenDetector.post_install import run
        run()
    except Exception as exc:
        print(
            f"\n[refgenDetector] WARNING: post-install setup failed: {exc}\n"
            "You can run it manually later with:\n"
            "    python -c 'from refgenDetector.post_install import run; run()'\n"
        )

class PostInstallCommand(install):
    def run(self):
        super().run()
        run_post_install()

class PostDevelopCommand(develop):
    def run(self):
        super().run()
        run_post_install()


# ── Setup ─────────────────────────────────────────────────────────────────────

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
            'refgenDetector-manager=refgenDetector.ref_manager:main'
        ],
    },
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    package_data={
        "refgenDetector": ["post_install.py"],   # ships post_install.py in the wheel
    }
)
