from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("ceq2_cy.pyx"))
