import io
import os
import sys
from distutils import sysconfig

import numpy as np
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools import setup, Extension

ROOT = os.path.dirname(__file__)


def get_version():
    # Borrowed this method from pysam
    sys.path.insert(0, "segmentation")
    # noinspection PyUnresolvedReferences,PyPackageRequirements
    import version
    return version.__version__


with io.open(os.path.join(ROOT, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

includes = [np.get_include(), sysconfig.get_config_var('INCLUDEDIR'), os.path.join(ROOT, 'src')]

CFLAGS = ['-Ofast']

setup(
    name="segmentation",
    packages=["segmentation"],
    version=get_version(),
    description='Mixture model segmentation',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Matej Usaj',
    author_email='m.usaj@utoronto.ca',
    url='https://github.com/usajusaj/segmentation',
    download_url='https://github.com/usajusaj/segmentation/archive/master.zip',
    keywords=['mixture', 'model', 'segmentation'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',

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
        'scipy',
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
                "segmentation/segmentation.pyx",
                "src/fastmm.cpp",
            ],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],  # suppress annoying compile warnings
            include_dirs=includes,
            extra_compile_args=CFLAGS,
            extra_link_args=CFLAGS
        )
    ])
)
