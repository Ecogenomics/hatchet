import os

from setuptools import setup, find_packages


def read_version():
    version_path = os.path.join(os.path.dirname(__file__), 'hatchet', 'VERSION')
    with open(version_path) as f:
        return f.readline().strip()


def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    with open(readme_path) as f:
        return f.read()


setup(
    name='hatchet',
    version=read_version(),
    author='Pierre-Alain Chaumeil',
    author_email='uqpchaum@uq.edu.au',
    maintainer='Pierre-Alain Chaumeil',
    description='Tools used to split the GTDB-Tk reference tree into smaller sub-trees.',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/Ecogenomics/hatchet',
    license='GPL3',
    packages=find_packages(exclude=['scripts', 'scripts.*']),
    package_data={'hatchet': ['VERSION']},
    scripts=['bin/hatchet'],
    python_requires='>=3.6',
    install_requires=[
        'dendropy',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)