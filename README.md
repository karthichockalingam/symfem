![Symfem](https://raw.githubusercontent.com/mscroggs/symfem/main/logo/logo.png)

[![Documentation status](https://readthedocs.org/projects/symfem/badge/?version=latest)](https://symfem.readthedocs.io/en/latest/?badge=latest)
[![Style checks](https://github.com/mscroggs/symfem/actions/workflows/style-checks.yml/badge.svg)](https://github.com/mscroggs/symfem/actions)
[![Run tests](https://github.com/mscroggs/symfem/actions/workflows/run-tests.yml/badge.svg)](https://github.com/mscroggs/symfem/actions)
[![Coverage Status](https://coveralls.io/repos/github/mscroggs/symfem/badge.svg?branch=main)](https://coveralls.io/github/mscroggs/symfem?branch=main)
[![PyPI](https://img.shields.io/pypi/v/symfem?color=blue&label=PyPI&logo=pypi&logoColor=white)](https://pypi.org/project/symfem/)
[![conda](https://img.shields.io/badge/Anaconda.org-2021.7.6-blue.svg?style=flat-square&label=conda&logo=anaconda)](https://anaconda.org/conda-forge/symfem)

Symfem is a symbolic finite element definition library, that can be used to
symbolically evaluate the basis functions of a finite element space.

## Installing Symfem
### Installing from repo
Symfem can be installed by downloading the [GitHub repo](https://github.com/mscroggs/symfem)
and running:

```bash
python3 setup.py install
```

### Installing using pip
The latest release of Symfem can be installed by running:

```bash
pip3 install symfem
```

### Installing using conda
The latest release of Symfem can be installed by running:

```bash
conda install symfem
```

## Testing Symfem
To run the Symfem unit tests, clone the repository and run:

```bash
python3 -m pytest test/
```

## Using Symfem
Documentation of the latest release version of Symfem can be found on
[Read the Docs](https://symfem.readthedocs.io/en/latest/).

## Contributing to Symfem
You can find information about how to contribute to Symfem [here](CONTRIBUTING.md).
