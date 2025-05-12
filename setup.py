from setuptools import setup, find_packages

#VERSION = '0.0.1' 
DESCRIPTION = 'SCOOTI'
LONG_DESCRIPTION = 'SCOOTI: Single Cell Optimization OBjective and Tradeoff Inference'



setup(
    name="SCOOTI",
    use_scm_version={"write_to": "SCOOTI/_version.py"},# Enables setuptools_scm to derive the version from SCM metadata
    setup_requires=["setuptools_scm"],  # Ensures setuptools_scm is installed during setup
    author="Chandrasekaran Lab, University of Michigan",
    author_email="csriram@umich.edu",
    description="SCOOTI: Single-Cell Objective Optimization and Tradeoff Inference",
    packages=find_packages(),
    install_requires=[
        "numpy==1.23.5",
        "pandas==1.5.3",
        "scikit-learn==1.1.3",
        "numba==0.56.4",
        "tqdm",
        "cobra==0.26.3",
        "scanpy==1.9.5",
        "seaborn==0.12.2",
        "phate==1.0.11",
        "adjustText==0.8",
        "torch==2.0.1"
    ],
    python_requires=">=3.8, <3.12",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # change if you use other
        "Operating System :: OS Independent",
    ],
)


"""
SCOOTI/
├── SCOOTI/
│   ├── __init__.py
│   ├── every_function.py
│   └── main.py
├── tests/
│   └── test_every_function.py
├── setup.py
├── README.md
├── LICENSE
└── MANIFEST.in   # (if you need to include non-code files)
"""
