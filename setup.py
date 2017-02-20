import os
import glob

tRNASeqTools_version='1.0'

requirements = [req.strip() for req in open('requirements.txt', 'rU').readlines() if not req.startswith('#')]

from setuptools import setup, find_packages 

if os.environ.get('USER','') == 'vagrant':
    del os.link

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

with open(os.path.join(os.path.dirname(__file__), 'README.md')) as readme:
    README = readme.read()

setup(
    name = "tRNASeqTools",
    version = tRNASeqTools_version,

    scripts = [script for script in glob.glob('bin/*')],
    include_package_data = True,

    packages = find_packages(),

    install_requires = requirements,

    author = "Meren Lab",
    author_email = "a.murat.eren@gmail.com",
    description = "Set of tools for tRNA-seq data analysis",
    license = "GPLv3+",
    keywords = "tRNA tRNA-seq bioinformatics",
    url = "https://merenlab.org/",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ],
)
