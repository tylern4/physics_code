from distutils.core import setup
from distutils.extension import Extension

import numpy
from Cython.Distutils import build_ext

include_dirs = [numpy.get_include(), "/usr/local/root/include", "../src"]
library_dirs = ["/usr/local/root/lib"]
extra_compile_args = [
    "-Wno-deprecated-register",
    "-Wno-#warnings",
    "-Wno-unused-variable",
    "-Wno-unused-function",
    "-Wno-sign-compare",
    "-pthread",
    "-std=c++14",
    "-m64"]

extra_link_args = ["-pthread",
                   "-m64"]

ext_modules = [
    # Extension("h10", ["h10.pyx"],
    #          language="c++",
    #          include_dirs=include_dirs,
    #          library_dirs=library_dirs,
    #          libraries=["Tree"],
    #          extra_compile_args=extra_compile_args,
    #          extra_link_args=extra_link_args),
    Extension("handeler", ["handeler.pyx"],
              language="c++",
              include_dirs=include_dirs,
              library_dirs=library_dirs,
              extra_compile_args=extra_compile_args,
              extra_link_args=extra_link_args)

]
setup(name='handeler',
      cmdclass={'build_ext': build_ext},
      ext_modules=ext_modules)
