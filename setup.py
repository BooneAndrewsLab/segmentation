import os
import sysconfig

import numpy as np
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools import setup, Extension

MKL = False

includes = [np.get_include(), sysconfig.get_config_var('INCLUDEDIR')]
for include_path in includes:
    if MKL:
        break

    for p, d, f in os.walk(include_path):
        if 'mkl.h' in f:
            MKL = True
            break

CFLAGS = ['-Ofast', '-std=c++17', '-Wcpp', '-DUSE_MKL']

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
            ],
            include_dirs=includes,
            extra_compile_args=CFLAGS,
            extra_link_args=CFLAGS
        )
    ])
)
