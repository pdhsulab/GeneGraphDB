from setuptools import find_packages, setup

setup(
    name="ggdbfetch",
    version="0.0.0_dev",
    description="A command line to retrieve contigs from GGDB.",
    url="https://github.com/pdhsulab/GeneGraphDB/GGDBFetch",
    author="Matt Durrant",
    author_email="mdurrant@berkeley.edu",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["click==8.0.1", "pyfastx", "pandas" "tqdm"],
    zip_safe=False,
    entry_points={"console_scripts": ["ggdbf = ggdbfetch.main:cli"]},
)
