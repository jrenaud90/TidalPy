# Contributing to TidalPy

Thank you for your interest in contributing to TidalPy! We welcome contributions from the community.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Enhancements](#suggesting-enhancements)

## Code of Conduct

This project adheres to a [code of conduct](https://tidalpy.readthedocs.io/en/latest/Overview/CoC.html). By participating, you are expected to uphold this code. Please be respectful and constructive in all interactions.

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check the existing issues to avoid duplicates. When you create a bug report, include as many details as possible:

- **Use a clear and descriptive title**
- **Describe the exact steps to reproduce the problem**
- **Provide specific examples** (code snippets, sample data)
- **Describe the behavior you observed** and what you expected
- **Include your environment details** (OS, Python version, TidalPy version)

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion:

- **Use a clear and descriptive title**
- **Provide a detailed description** of the suggested enhancement
- **Explain why this enhancement would be useful**
- **Include examples** of how the feature would be used

### Pull Requests

Please feel free to make pull requests! Also don't hesitate to make draft PRs so the developer can assist with your changes.

1. Fork the repo and create your branch from `main`
2. If you've added code, please add tests
3. If you've changed APIs, please update the documentation
4. Ensure the test suite passes by running `pytest Tests\` (recommended you use `pytest-xdist` and use multiple threads; there are a lot of tests!).
5. Submit your pull request!

## Getting Started

### Prerequisites

- Python >= 3.9
- A C and C++ compiler that supports C++20 (MSVC on Windows, GCC on Linux, Apple's clang on MacOS)
- Git
- A GitHub account

### Setting Up Your Development Environment

1. **Fork and clone the repository** (with its submodules, the header-only C++ libraries in `Dependencies/`):
   ```bash
   git clone --recursive https://github.com/YOUR-USERNAME/TidalPy.git
   cd TidalPy
   ```
   For a clone made without `--recursive`, run `git submodule update --init`.

2. **Create a virtual environment:**
   ```bash
   python -m venv venv # Or conda environment if you prefer
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Build and install TidalPy with its development dependencies:**
   ```bash
   pip install -v ".[dev]"
   ```
   Reinstall after every change: a change to a `.pyx`, `.pxd`, or `.hpp` file needs the extensions recompiled, which the reinstall does.

4. **Create a new branch for your feature:**
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Workflow

1. **Make your changes** in your feature branch
2. **Write or update tests** for your changes
3. **Run the test suite** to ensure everything passes
4. **Update documentation** if needed
5. **Commit your changes** with clear, descriptive messages
6. **Push to your fork** and submit a pull request

### Commit Message Guidelines

- Use the present tense ("Add feature" not "Added feature")
- Limit the first line to 72 characters or less
- Reference issues and pull requests liberally after the first line

Example:
```
Add rheological model for Maxwell material

- Implement frequency-dependent compliance
- Add tests for various temperature ranges
- Update documentation with usage examples

Fixes #123
```

## Repository Layout

- `TidalPy/<Module>/`: the source of each module (`Structures`, `Material`, `Rheology`, `Viscosity`, `PartialMelt`, `Cooling`, `Radiogenics`, `RadialSolver`, `Tides`, `Dynamics`, `Stellar`, `Utilities`). C++ headers (`*_.hpp`) hold the physics; Cython files (`*.pyx`, `*.pxd`) wrap them; Python holds configuration, file I/O, and plotting. `TidalPy/WorldPack/` holds the bundled world and system TOML files, and `TidalPy/defaultc.py` the default configuration.
- `Tests/Test_<Module>/`: the tests of each module, plus `Tests/Test_E2E/` (configuration-to-physics runs) and `Tests/Test_Package/` (import, configuration, and logging).
- `Documentation/`: the Sphinx (myst-parser) pages, one folder per module.
- `Demos/` and `Benchmarks/`: tutorial notebooks and validation or performance notebooks.
- `Dependencies/`: git submodules for the header-only C++ libraries (Eigen, xsf, spdlog).
- `setup.py` and `cython_extensions.json`: the build script and the list of compiled extensions.

## Coding Standards

### Python Style Guide

- Follow [PEP 8](https://pep8.org/) style guidelines.
- Follow [Numpy](https://numpydoc.readthedocs.io/en/latest/format.html) style for function and class docstrings.
- Use meaningful variable and function names.
- Try to keep to a maximum line length of 120 characters.
- Use type hints where appropriate.

### Code Quality Tools

We use the following tools to maintain code quality:

- **Linting:** `ruff`

Run these before submitting:
```bash
ruff check .
```

### Documentation Style

- Use NumPy-style docstrings
- Include parameter types and descriptions
- Provide usage examples for public APIs
- Keep documentation up-to-date with code changes

Example:
```python
def calc_tidal_heating(world, orbital_frequency, eccentricity, host_mass):
    """
    Calculate the tidal heating of a world on a synchronous orbit.

    Parameters
    ----------
    world : BaseWorld
        The world, with its equation of state solved.
    orbital_frequency : float
        Orbital mean motion [rad s-1].
    eccentricity : float
        Orbital eccentricity.
    host_mass : float
        Mass of the tidal host [kg].

    Returns
    -------
    float
        Tidal heating [W].

    Examples
    --------
    >>> heating = calc_tidal_heating(europa, orbital_frequency, 0.009, mass_jupiter)
    >>> print(f"Heating: {heating:.2e} W")
    """
```

## Testing

### Running Tests

TidalPy has lots of tests! It is highly recommended you install `pip install pytest-xdist` and use multiple cores with `pytest -n logical Tests/`. Prefer `-n logical` over `-n auto`. When `psutil` is installed, `auto` counts only physical cores, which is half the workers on a machine with hyper-threading, `-n logical` overcomes this limitation.

The tests import the installed TidalPy, not the source tree: the repository's `conftest.py` removes the repository root from the import path. Run `pytest` from the repository root after reinstalling. `Tests/conftest.py` points `TIDALPY_DATA_DIR` at a fresh temporary directory before TidalPy is imported, so the suite always runs with the packaged configuration, worlds, and materials, never reads or writes your own data directory, and gives the same results whatever your `TidalPy_Configs.toml` says.

```bash
# Run all tests
pytest Tests/  # Or with the added `-n logical` flag.

# Run specific test file
pytest Tests/Test_Rheology/test_rheology_01.py
```

_Note that multiple warnings may show while you are running tests. These are likely normal warnings and are expected. TidalPy will try to raise an `Exception` (which `pytest` should catch automatically) when there is a serious problem._

### Writing Tests

- Write unit tests for new functions and classes
- Include edge cases and error conditions
- Use descriptive test names that explain what is being tested
- Place tests in the `Tests/Test_<Module>/` directory of the module they test

Example:
```python
from TidalPy.Rheology import Maxwell


def test_maxwell_zero_frequency():
    """The Maxwell complex modulus vanishes at zero forcing frequency."""
    model = Maxwell()
    complex_modulus = model.calc_complex_modulus(5.0e10, 1.0e21, 0.0)  # (modulus [Pa], viscosity [Pa s], frequency [rad s-1])
    assert abs(complex_modulus) < 1.0e-6 * 5.0e10
```

## Documentation

Documentation is built using Sphinx and hosted at [https://tidalpy.readthedocs.io/en/latest/](https://tidalpy.readthedocs.io/en/latest/).

### Adding Documentation

- Update the relevant `.md` pages in the `Documentation/<Module>/` folder, and the toctree in that folder's `index.md` when adding a page
- Run every code example you add or change against an installed build, and check that every relative link resolves
- Include docstrings in your code
- Add examples and tutorials for new features; notebooks in `Demos/` and `Benchmarks/` are copied into the documentation build and shown with their stored outputs

## Submitting Changes

### Pull Request Process

1. **Update the CHANGES.md** with details of your changes
2. **Ensure all tests pass** and coverage doesn't decrease
3. **Update documentation** to reflect any changes
4. **Request review** from maintainers
5. **Address feedback** from reviewers
6. Once approved, a maintainer will merge your PR

### Pull Request Checklist

- [ ] Code follows the project's style guidelines
- [ ] Tests added/updated and all tests pass
- [ ] Documentation updated
- [ ] CHANGELOG updated
- [ ] Commit messages are clear and descriptive
- [ ] Branch is up-to-date with main

### Marking Issues or ToDos Inside the Source Code
Ideally you will not have to leave todos or mark issues in the code. But if you do please follow the format below:

- "TODO: message" - For items that should be done in the future, perhaps there is a question about implementation
  details.
- "FIXME: message" - Critical issues that should be fixed immediately. These should only be used during PR drafting
  and should never make it into a release.
- "OPT: message" - Areas of code you think can be optimized but are deferring it for the future.

**All of these should be turned into GitHub issues! Marking them in the code should only serve as an additional marker.**

## Questions?

If you have questions about contributing, feel free to:

- Open an issue with the `question` label
- Email the team: [TidalPy@gmail.com](mailto:TidalPy@gmail.com)
    - Also feel free to ask us to invite you to TidalPy's Slack account for faster communications.
- Check the documentation at [TidalPy.info](https://TidalPy.info)

## License

By contributing to TidalPy, you agree that your contributions will be licensed under the CC-BY-SA-4.0 License.

---

Thank you for contributing to TidalPy!
