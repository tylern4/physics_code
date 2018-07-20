from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(ext_modules=[Extension("physics",
                             ["physics.pyx"],
                             language="c++",
                             include_dirs=["/usr/local/root/include"],
                             library_dirs=["/usr/local/root/lib"],
                             libraries=["Tree"],
                             extra_compile_args=[
                                 "-Wno-unused-variable", "-Wno-unused-function",
                                 "-Wno-sign-compare", "-D__LZ4__", "-pthread",
                                 "-stdlib=libc++", "-std=c++1z", "-m64"],
                             extra_link_args=["-D__LZ4__", "-pthread",
                                              "-stdlib=libc++", "-std=c++1z", "-m64"])],
      cmdclass={'build_ext': build_ext})
