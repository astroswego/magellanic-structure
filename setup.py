#!/usr/bin/env python3
"""magellanic-structure: description goes here

long description goes here
"""

DOCLINES = __doc__.split('\n')

CLASSIFIERS = """\
Programming Language :: Python
Programming Language :: Python :: 3
Intended Audience :: Science/Research
"""

MAJOR = 0
MINOR = 1
MICRO = 0
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

def setup_package():
    metadata = dict(
        name = 'magellanic-structure',
        url = 'https://github.com/astroswego/magellanic-structure',
        description = DOCLINES[0],
        long_description = "\n".join(DOCLINES[2:]),
        version = VERSION,
        package_dir = {'': 'src'},
        packages = [
            'magstruct',
            'magstruct_scripts'
        ],
        entry_points = {
            'console_scripts': [
                'magstruct = magstruct_scripts.magstruct:main'
            ]
        },
        classifiers = [f for f in CLASSIFIERS.split('\n') if f],
        requires = [
            'numpy (>= 1.8.0)',
            'sklearn (>= 0.14.0)'
        ]
    )

    from setuptools import setup

    setup(**metadata)

if __name__ == '__main__':
    setup_package()
