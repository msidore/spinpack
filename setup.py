from setuptools import setup, Extension
import numpy.distutils.misc_util

setup(name='spinpack',
    packages=['spinpack'],
    ext_modules=[Extension("spinpack.spin", ["spinpack/spin.c"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    version='0.2',
    description='Get the spinangle from a NAMD trajectory',
    classifiers=[
        'Development Status :: 0 - Alpha',
        'License :: Friand_Dur',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    author='Marlon Sidore',
    author_email='marlon.sidore@gmail.com',
    license='CRAPL',
    install_requires=[
        'MDAnalysis',
        'numpy',
    ],
      zip_safe=False)
