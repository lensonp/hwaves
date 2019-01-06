from setuptools import setup, find_packages

longdesc = 'Current capabilities: '\
    'Compute and plot radial wavefunction density; '\
    'Compute and plot spherical harmonics; '\
    'Compute real-space density of hydrogenic eigenfunctions; '\
    'Compute real-space density of real-valued wavefunctions; '\
    'Plot (and provide polygon data) for real-space density isosurfaces of any value'

setup(
    name='hwaves',
    version='0.0.3',
    url='https://github.com/lensonp/hwaves.git',
    description='Compute and plot hydrogenic wavefunctions',
    long_description=longdesc,
    author='Lenson A. Pellouchoud',
    license='BSD',
    author_email='',
    install_requires=['numpy','scipy','masstable','periodictable'],
    packages=find_packages(),
    package_data={}
    )


