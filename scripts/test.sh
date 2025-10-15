#!/usr/bin/env bash
# Test runner script

set -e

echo "Running viewtools test suite..."

# Check if pytest is available
if ! command -v pytest &> /dev/null; then
    echo "pytest not found. Installing development dependencies..."
    pip install -e .[dev]
fi

echo ""
echo "=== Running unit tests ==="
python -m pytest tests/unit/ -v --cov=viewtools --cov-report=term-missing

echo ""
echo "=== Running integration tests ==="
python -m pytest tests/integration/ -v

echo ""
echo "=== Running all tests with coverage ==="
python -m pytest tests/ -v --cov=viewtools --cov-report=html --cov-report=term-missing

echo ""
echo "=== Test Summary ==="
echo "✓ Unit tests completed"
echo "✓ Integration tests completed" 
echo "✓ Coverage report generated in htmlcov/"

echo ""
echo "=== Code Quality Checks ==="

# Run black formatter check
if command -v black &> /dev/null; then
    echo "Checking code formatting with black..."
    black --check viewtools/ tests/ || echo "⚠️  Code formatting issues found. Run 'black viewtools/ tests/' to fix."
else
    echo "black not found, skipping format check"
fi

# Run ruff linter
if command -v ruff &> /dev/null; then
    echo "Running linter checks with ruff..."
    ruff check viewtools/ tests/ || echo "⚠️  Linting issues found. Run 'ruff check --fix viewtools/ tests/' to fix."
else
    echo "ruff not found, skipping lint check"
fi

echo ""
echo "Tests completed successfully! 🎉"
