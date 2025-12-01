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

This project adheres to a [code of conduct](https://tidalpy.readthedocs.io/en/latest/CodeOfConduct.html). By participating, you are expected to uphold this code. Please be respectful and constructive in all interactions.

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

- Python >= 3.8
- Git
- A GitHub account

### Setting Up Your Development Environment

1. **Fork and clone the repository:**
   ```bash
   git clone https://github.com/YOUR-USERNAME/TidalPy.git
   cd TidalPy
   ```

2. **Create a virtual environment:**
   ```bash
   python -m venv venv # Or conda environment if you prefer
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install development dependencies:**
   ```bash
   pip install -e ".[dev]"
   ```

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
def calculate_tidal_heating(body, orbit, rheology = Maxwell):
    """
    Calculate tidal heating for a planetary body.

    Parameters
    ----------
    body : Body
        The planetary body object
    orbit : Orbit
        The orbital parameters
    rheology : Rheology, default = Maxwell
        The rheological model to use

    Returns
    -------
    float
        Tidal heating in Watts

    Examples
    --------
    >>> heating = calculate_tidal_heating(europa, orbit, maxwell)
    >>> print(f"Heating: {heating:.2e} W")
    """
```

## Testing

### Running Tests

TidalPy has lots of tests! It is highly recommended you install `pip install pytest-xdist` and use multiple cores
with `pytest -n auto Tests/`

```bash
# Run all tests
pytest Tests/

# Run specific test file
pytest Tests/test_rheology.py
```

### Writing Tests

- Write unit tests for new functions and classes
- Include edge cases and error conditions
- Use descriptive test names that explain what is being tested
- Place tests in the `Tests/` directory

Example:
```python
def test_maxwell_rheology_zero_frequency():
    """Test that Maxwell rheology returns correct value at zero frequency."""
    model = MaxwellRheology(viscosity=1e21, shear_modulus=5e10)
    result = model.compliance(frequency=0)
    assert result > 0
```

## Documentation

Documentation is built using Sphinx and hosted at [https://tidalpy.readthedocs.io/en/latest/](https://tidalpy.readthedocs.io/en/latest/).

### Adding Documentation

- Update relevant `.rst` or `.md` files in the `Documentation/` directory
- Include docstrings in your code
- Add examples and tutorials for new features

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

## Questions?

If you have questions about contributing, feel free to:

- Open an issue with the `question` label
- Email the maintainer: joe.p.renaud@gmail.com
- Check the documentation at [TidalPy.info](https://TidalPy.info)

## License

By contributing to TidalPy, you agree that your contributions will be licensed under the CC-BY-SA-4.0 License.

---

Thank you for contributing to TidalPy!