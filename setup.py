from setuptools import setup, find_packages

setup(
    name="genegraphdb",
    version='0.0.0_dev',
    description='A command line tool to create and query a gene graph databases.',
    url='https://github.com/pdhsulab/GeneGraphDB',
    author="Matt Durrant",
    author_email="mdurrant@berkeley.edu",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click==8.0.1',
        'neo4j==4.3.1',
        'pytz==2021.1',
        'biopython==1.79',
        'numpy==1.21.0'
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'ggdb = genegraphdb.main:cli'
        ]
}
)