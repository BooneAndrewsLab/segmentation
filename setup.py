import sysconfig

import numpy as np
from Cython.Distutils import build_ext
from setuptools import setup, Extension

our_flags = ["-march=native", "-std=c++11", "-DNDEBUG", "-Ofast", "-Wcpp"]

ext_modules = [
    Extension("segmentation.segmentation", ["segmentation/segmentation.pyx", ],
              language="c++",
              include_dirs=[np.get_include(), sysconfig.get_config_var('INCLUDEDIR')],
              extra_compile_args=our_flags,
              extra_link_args=["-std=c++11"],
              ),
]

setup(
    name="segmentation",
    packages=["segmentation"],
    version='0.1.2',
    description='Mixture model segmentation',
    author='Matej Usaj',
    author_email='m.usaj@utoronto.ca',
    url='https://github.com/usajusaj/segmentation',
    download_url='https://github.com/usajusaj/segmentation/archive/master.zip',
    keywords=['mixture', 'model', 'segmentation'],
    classifiers=[],
    cmdclass={"build_ext": build_ext},
    ext_modules=ext_modules,
    install_requires=['cython', 'scikit-image', 'numpy'],
)
