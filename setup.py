import sysconfig

import numpy as np
from Cython.Distutils import build_ext
from setuptools import setup, Extension

our_flags = ["-march=native", "-std=c++11", "-DNDEBUG", "-Ofast", "-Wcpp"]
extra_compile_args = sysconfig.get_config_var('CFLAGS').split()

for f in extra_compile_args[:]:
    if f.startswith('-march=') or f.startswith('-std=') or f.startswith('-O'):
        extra_compile_args.remove(f)

extra_compile_args += our_flags

ext_modules = [
    Extension("segmentation.segmentation", ["segmentation/segmentation.pyx", ],
              language="c++",
              include_dirs=[np.get_include()],
              extra_compile_args=extra_compile_args,
              extra_link_args=["-std=c++11"],
              ),
]

setup(
    name="segmentation",
    cmdclass={"build_ext": build_ext},
    ext_modules=ext_modules)
