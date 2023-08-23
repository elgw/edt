# python3 eudist_setup.py build_ext --inplace
# python3 eudist_setup.py install

from distutils.core import setup, Extension

module1 = Extension('eudist',
                    sources = ['eudist_py.c'])

setup (name = 'eudist',
       version = '1.0',
       description = 'Euclidean distance transform',
ext_modules = [module1])
