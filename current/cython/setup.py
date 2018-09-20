from distutils.core import setup
from distutils.extension import Extension

import numpy
from Cython.Distutils import build_ext

include_dirs = [numpy.get_include()] + ["/usr/local/root/include"]

setup(ext_modules=[Extension("h10",
                             ["h10.pyx"],
                             language="c++",
                             include_dirs=include_dirs,
                             library_dirs=["/usr/local/root/lib"],
                             libraries=["Tree"],
                             extra_compile_args=[
                                 "-Wno-deprecated-register",
                                 "-Wno-#warnings",
                                 "-Wno-unused-variable",
                                 "-Wno-unused-function",
                                 "-Wno-sign-compare",
                                 "-pthread",
                                 "-stdlib=libc++",
                                 "-std=c++14",
                                 "-m64"],
                             extra_link_args=["-pthread",
                                              "-stdlib=libc++",
                                              "-std=c++14",
                                              "-m64"])],
      cmdclass={'build_ext': build_ext})
