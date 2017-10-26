from distutils.core import setup, Extension
from Cython.Build import cythonize
import os
import pkg_resources

pysam = pkg_resources.resource_filename('pysam', '')

# https://docs.python.org/2/distutils/apiref.html
pysamx = Extension(name="pysamx",
                   sources=["pysamx.pyx"],
                   include_dirs=[os.path.join(pysam, 'include/htslib'), pysam],
                   library_dirs=[pysam],
                   libraries=["calignmentfile", "chtslib", "calignedsegment"],
                   extra_compile_args=["-I {} -I {}".format(os.path.join(pysam, 'include/htslib'), pysam)])

setup(
    ext_modules = cythonize([pysamx])
)
