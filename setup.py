import setuptools

setuptools.setup(
    name="IFAA",
    version="0.1.0",
    url="https://github.com/lovestat/ssfunc",
    author="Shangchen Song",
    author_email="s.song@ufl.edu",
    description="Allows conversion of Roman numerals to ints (and vice versa)",
    long_description=open('README.md').read(),
    packages=setuptools.find_packages(),
    install_requires=[],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
    include_package_data=True,
    package_data={'': ['data/*.csv']},
)
