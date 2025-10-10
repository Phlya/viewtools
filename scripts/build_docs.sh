#!/bin/bash
# Build documentation with Sphinx

set -e

echo "Building documentation..."

# Check if we're in the right directory
if [ ! -f "pyproject.toml" ]; then
    echo "Error: Run this script from the project root directory"
    exit 1
fi

# Install development dependencies if needed
echo "Checking dependencies..."
pip install -e .[dev] > /dev/null 2>&1 || {
    echo "Failed to install dependencies"
    exit 1
}

# Create docs directory structure if needed
mkdir -p docs/_static
mkdir -p docs/_templates

# Change to docs directory
cd docs/

# Clean previous builds
echo "Cleaning previous builds..."
rm -rf _build/

# Build HTML documentation
echo "Building HTML documentation..."
sphinx-build -b html . _build/html

echo ""
echo "Documentation built successfully!"
echo "Open docs/_build/html/index.html to view the documentation."

# Check if we can serve it
if command -v python3 &> /dev/null; then
    echo ""
    echo "To serve documentation locally, run:"
    echo "  cd docs/_build/html && python3 -m http.server 8000"
fi
