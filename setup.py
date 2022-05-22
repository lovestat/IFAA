import setuptools

setuptools.setup(
    name="IFAA",
    version="0.1.0",
    url="https://github.com/lovestat/IFAA",
    author="Zhigang Li, Quran Wu, Shangchen Song",
    author_email="lzg2151@gmail.com, s.song@ufl.edu",
    description="Robust association identification and inference for absolute abundance in microbiome analyses",
    long_description=open('README.md').read(),
    packages=setuptools.find_packages(),
    install_requires=['joblib>=0.10.3', 'scipy', 'numpy', 'pandas', 'glmnet_py', 'matplotlib', 
                      'tqdm', 'statsmodels'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    include_package_data=True,
    package_data={'': ['data/*.csv']},
)
