from setuptools import setup, find_packages

requires = [
    'click',
    'biopython>=1.7',
    'pandas'
    ]

setup(
    name="mbarq",
    version="1.0",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description=("Bioinformatic toolkit for barcode mapping and counting"),
    license="LICENSE",
    keywords="mbarq",
    install_requires=requires,
    #url = "https://github.com/SushiLab/TNSEQ_DEV",
    packages=find_packages(),
    entry_points={
        'console_scripts': ['mbarq=mbarq.cli:main'],
    }
)
