import os
from distutils import sysconfig

import numpy as np
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools import setup, Extension

ROOT = os.path.dirname(os.path.abspath(__file__))

includes = [np.get_include(), sysconfig.get_config_var('INCLUDEDIR'), os.path.join(ROOT, 'src')]

CFLAGS = ['-Ofast']

setup(
    name="segmentation",
    packages=["segmentation"],
    version='0.1.3',
    description='Mixture model segmentation',
    author='Matej Usaj',
    author_email='m.usaj@utoronto.ca',
    url='https://github.com/usajusaj/segmentation',
    download_url='https://github.com/usajusaj/segmentation/archive/master.zip',
    keywords=['mixture', 'model', 'segmentation'],
    classifiers=[
        'Development Status :: 3 - Beta',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    setup_requires=[
        'cython',
        'numpy'
    ],
    install_requires=[
        'numpy',
        'scikit-image'
    ],
    entry_points={
        'console_scripts': [
            'segment=segmentation.cli:main',
        ],
    },
    cmdclass={"build_ext": build_ext},
    ext_modules=cythonize([
        Extension(
            name="segmentation.segmentation",
            language="c++",
            sources=[
                "src/segmentation.pyx",
                "src/fastmm.cpp",
            ],
            define_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],  # suppress annoying compile warnings
            include_dirs=includes,
            extra_compile_args=CFLAGS,
            extra_link_args=CFLAGS
        )
    ])
)
