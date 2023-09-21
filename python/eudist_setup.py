# python3 eudist_setup.py build_ext --inplace
# python3 eudist_setup.py install

import sysconfig

from distutils.core import setup, Extension

extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
extra_compile_args += ["-DNDEBUG", "-O3"]

module1 = Extension('eudist.eudist',
                    include_dir=['../src/'],
                    sources = ['../src/eudist_py.c', '../src/eudist.c'],
                    extra_compile_args=extra_compile_args)

setup (name = 'eudist',
       version = '0.0.1',
       description = 'Euclidean distance transform',
ext_modules = [module1])
