# Contributing

## Development Setup

1. Clone the repository:

```bash
git clone https://github.com/phlya/viewtools.git
cd viewtools
```

2. Install in development mode:

```bash
pip install -e .[dev]
```

3. Run tests:

```bash
pytest
```

## Code Style

We use:

- `black` for code formatting
- `ruff` for linting
- Type hints for all functions

Format code:

```bash
black viewtools/
ruff check viewtools/
```

## Testing

Run the test suite:

```bash
pytest tests/
```

With coverage:

```bash
pytest --cov=viewtools tests/
```

## Documentation

Build documentation locally:

```bash
cd docs/
sphinx-build -b html . _build/html
```

Or use the build script:

```bash
python scripts/build_docs.py
```

The documentation will be built in `docs/_build/html/`.

## Pull Requests

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if needed
5. Update documentation
6. Submit a pull request

## Release Process

1. Update version in `pyproject.toml`
2. Update `CHANGELOG.md`
3. Create a git tag
4. Push to GitHub
5. Create a release on GitHub
