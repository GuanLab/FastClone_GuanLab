#!/usr/bin/env python3

import setuptools

with open('README.md', 'r') as fh:
        long_description = fh.read()

setuptools.setup(
    name='fastclone-guanlab',
    version='1.0.8',
    description='An inference tool for tumour subclonal composition',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Hongjiu Zhang, Yuanfang Guan, Yao Xiao',
    author_email='zhanghj@umich.edu, gyuanfan@umich.edu, yxiaoy@umich.edu',
    url='http://github.com/GuanLab/FastClone_GuanLab',
    license='GPLv3+',
    keywords='bioinformatics',
    python_requires='>=3.5.0',
    packages=setuptools.find_packages(exclude=['tests']),
    entry_points={
        'console_scripts': [
            'fastclone = fastclone.__main__:main',
        ],
    },
    package_data={
        'fastclone': ['test/data/*'],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        ('License :: OSI Approved '
         ':: GNU General Public License v3 or later (GPLv3+)'),
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    #setup_requires=['pytest-runner'],
    #tests_require=['pytest', 'pytest-datadir'],
    install_requires=['logbook', 'pandas', 'scikit-learn', 'fire', 'networkx', 'matplotlib'],
)
