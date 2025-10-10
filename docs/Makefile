# Minimal makefile for Sphinx documentation

SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx-build
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

# Specific targets
html:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS)

clean:
	rm -rf $(BUILDDIR)/*

livehtml:
	@if command -v sphinx-autobuild >/dev/null 2>&1; then \
		sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html"; \
	else \
		echo "sphinx-autobuild not found. Install with: pip install sphinx-autobuild"; \
		exit 1; \
	fi

.PHONY: html clean livehtml
