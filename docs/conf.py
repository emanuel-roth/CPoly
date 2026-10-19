# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
try:
    import sage.all
except ImportError:
    sage = None

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + "/.."))

# -- Project information -----------------------------------------------------

project = "CPoly"
copyright = "2026"
author = "Emanuel Roth"
# The full version, including alpha/beta/rc tags
release = "1.0"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinxcontrib.inkscapeconverter",
]

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"
mathjax3_config = {
    "tex": {
        "inlineMath": [["$", "$"], ["\\(", "\\)"]],
        "displayMath": [["$$", "$$"], ["\\[", "\\]"]],
        "linebreaks": {
            "automatic": True
        }
    }
}

latex_engine = "lualatex"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_title = "CPoly"
html_theme = "furo"
html_static_path = ["_static"]
html_css_files = [
    "custom.css",
]
html_theme_options = {
    "light_logo": "cpolylogo.svg",
    "dark_logo": "cpolylogodark.svg",
}
html_favicon = "_static/favicon.ico"
html_show_sourcelink = False
viewcode_follow_imported_members = False
autodoc_default_options = {
    "members": True,
    "inherited-members": True,
    "show-inheritance": True,
}