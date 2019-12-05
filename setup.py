#!/usr/bin/env python
from setuptools import setup, find_packages

install_requires = """
    astropy<3.0.0
    numpy>=1.3.0,<1.17.0
    psycopg2
    python-casacore
    python-dateutil>=1.4.1
    pytz
    scipy>=0.7.0,<1.3.0
    sqlalchemy>=1.0.0
    alembic
    monotonic
    astroML
    matplotlib<3.0
    pandas
    astroquery
    tkp
    """.split()

# extras_require = {
#     'monetdb': ['sqlalchemy_monetdb>=0.9.1'],
# }

tkp_scripts = [
    "traptools/bin/GetThumbnails.py",
    ]

# package_data = {'tkp': [
#     'config/*/*',
#     'db/sql/statements/batch',
#     'db/sql/statements/*/*.sql'
# ]}

package_list = find_packages(where='.', exclude=['tests'])

setup(
    name="traptools",
    version="0.1",
    packages=package_list,
    scripts=tkp_scripts,
    # package_data=package_data,
    description="TraP Tools",
    author="TKP",
    author_email="discovery@transientskp.org",
    url="http://docs.transientskp.org/",
    install_requires=install_requires,
    # extras_require=extras_require
)