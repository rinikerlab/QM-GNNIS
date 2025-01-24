from setuptools import setup, find_packages

setup(
    name="implicitml", 
    version="0.1.0",  
    description="A package for implicit machine learning models",
    author="Paul Katzberger, Felix Pultar",  
    author_email="paul.katzberger@phys.chem.ethz.ch",  
    url="https://github.com/rinikerlab/implicitml",  
    packages=find_packages(), 
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # License of your project
        "Operating System :: OS Independent",
    ],
)