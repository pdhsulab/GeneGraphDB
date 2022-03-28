from setuptools import setup, find_packages

setup(
    name="annotategff",
    version='0.0.0_dev',
    description='A command line to annotate gffs for contigs from GGDB.',
    url='https://github.com/pdhsulab/GeneGraphDB/GGDBFetch',
    author="Jacob Luo",
    author_email="jluo5972@berkeley.edu",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click==8.0.1',
        'pandas',
        'biopython==1.79',
        'numpy==1.21.0',
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'ggdbf = ggdbfetch.main:cli'
        ]
}
)