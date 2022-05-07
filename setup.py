import setuptools

setuptools.setup(
    name="IFAA",
    version="0.1.0",
    url="https://github.com/lovestat/IFAA",
    author="Shangchen Song",
    author_email="s.song@ufl.edu",
    description="Robust association identification and inference for absolute abundance in microbiome analyses",
    long_description=open('README.md').read(),
    packages=setuptools.find_packages(),
    install_requires=['joblib>=0.10.3'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
    include_package_data=True,
    package_data={'': ['data/*.csv']},
)
