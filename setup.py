from distutils.core import setup, Extension
import glob
from pysam import get_include as pysam_get_include


try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    raise ImportError("Requires cython to "
            "be installed before running setup.py (pip install cython)")
try:
    import numpy as np
except ImportError:
    raise ImportError("Requires numpy to "
            "be installed before running setup.py (pip install numpy)")
try:
    import pysam
except ImportError:
    raise ImportError("Requires pysam to "
            "be installed before running setup.py (pip install pysam)")
try:
    import scipy
except ImportError:
    raise ImportError("Requires scipy to "
            "be installed before running setup.py (pip install scipy)")
try:
    import sklearn
except ImportError:
    raise ImportError("Requires sklearn to "
            "be installed before running setup.py (pip install sklearn)")
try:
    import pandas
except ImportError:
    raise ImportError("Requires pandas to "
            "be installed before running setup.py (pip install pandas)")


include_path = ['/usr/include']
include_path.extend(pysam_get_include())
include_path.append('/stor/work/Lambowitz/cdw2854/src/miniconda3/include')
ext_modules=cythonize([
        Extension('*', ['tgirt_smRNA_correction/*.pyx'],
                  include_dirs = list(include_path) )
])


setup(
    name='tgirt_smRNA_correction',
    version='0.1',
    description='Coverage correction for TGIRT-seq data',
    url='',
    author='Douglas Wu',
    author_email='wckdouglas@gmail.com',
    license='MIT',
    packages=['tgirt_smRNA_correction'],
    zip_safe=False,
    scripts = glob.glob('bin/*'),
    ext_modules = ext_modules,
    install_requires=[
          'cython',
          'numpy',
          'pysam>=0.11.0',
          'pandas>=0.20.2',
          'tqdm',
          'scipy>=0.19.0'
      ],
    package_data={'tgirt_smRNA_correction': ['model/*.pkl']},
    cmdclass = {'build_ext': build_ext}
)
