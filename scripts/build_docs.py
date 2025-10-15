#!/usr/bin/env python3
"""Build documentation using Sphinx."""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def run_command(cmd, cwd=None):
    """Run a command and handle errors."""
    try:
        result = subprocess.run(
            cmd, shell=True, cwd=cwd, check=True, capture_output=True, text=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        print(f"Exit code: {e.returncode}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)


def main():
    """Main build function."""
    # Get project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    docs_dir = project_root / "docs"

    print("Building viewtools documentation...")

    # Change to project root
    os.chdir(project_root)

    # Install dependencies
    print("Installing dependencies...")
    run_command("pip install -e .[dev]")

    # Create necessary directories
    (docs_dir / "_static").mkdir(exist_ok=True)
    (docs_dir / "_templates").mkdir(exist_ok=True)

    # Clean previous build
    build_dir = docs_dir / "_build"
    if build_dir.exists():
        print("Cleaning previous build...")
        shutil.rmtree(build_dir)

    # Build documentation
    print("Building HTML documentation...")
    run_command(f"sphinx-build -b html {docs_dir} {build_dir}/html")

    print("\nDocumentation built successfully!")
    print(f"Open {build_dir}/html/index.html to view the documentation.")

    # Offer to serve locally
    print("\nTo serve documentation locally, run:")
    print(f"  cd {build_dir}/html && python3 -m http.server 8000")


if __name__ == "__main__":
    main()
