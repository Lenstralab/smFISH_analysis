import setuptools
import os


with open('README.md', 'r') as fh:
    long_description = fh.read()


setuptools.setup(
    name='smfish',
    version='2021.12.0',
    author='Lenstra lab NKI',
    author_email='t.lenstra@nki.nl',
    description='Single molecule FISH code for the Lenstra lab.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/Lenstralab/smFISH_analysis',
    packages=['smfish'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
    tests_require=['pytest-xdist'],
    install_requires=['numpy', 'pandas', 'scipy', 'tqdm', 'tifffile', 'matplotlib', 'pyyaml', 'seaborn', 'scikit-image',
                      'parfor', 'lmfit', 'psutil', 'scikit-learn', 'rtree',
                      'tllab_common@git+https://github.com/Lenstralab/tllab_common.git@dbfe63427cabb3d9409468592edabe17a1dd27be'],
    scripts=[os.path.join('bin', script) for script in
             os.listdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin'))],
)
