# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import sphinx_rtd_theme
from ase.utils.sphinx import mol_role

sys.path.insert(os.path.abspath('../..'))
sys.path.insert(0, os.path.abspath("../hetbuilder"))
sys.path.insert(0, os.path.abspath(".."))

import re

VERSIONFILE = "../hetbuilder/__init__.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

import recommonmark
from recommonmark.transform import AutoStructify

# -- Project information -----------------------------------------------------

project = "hetbuilder"
copyright = "2021, Roman Kempt"
author = "Roman Kempt"

# The full version, including alpha/beta/rc tags
release = version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.


intersphinx_mapping = {
    "ase": ("https://wiki.fysik.dtu.dk/ase", None),
    "python": ("https://docs.python.org/3.7", None),
}

extensions = [
#    "nbsphinx",
    "sphinx.ext.autodoc",
#    "sphinx_markdown_tables",
    "sphinx.ext.todo",
    "sphinx.ext.githubpages",
    "recommonmark",
    "sphinxcontrib.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "misc", "**.ipynb_checkpoints"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
# html_static_path = ["_static"]

pygments_style = "sphinx"
master_doc = "index"

# -- Extension configuration -------------------------------------------------

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "restructuredtext",
    ".md": "markdown",
}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = False
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = False
# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True


# At the bottom of conf.py
github_doc_root = "https://github.com/rtfd/recommonmark/tree/master/doc/"

# dummv variables to avoid errors
nbsphinx_allow_errors = True


def setup(app):
    app.add_config_value(
        "recommonmark_config",
        {
            "url_resolver": lambda url: github_doc_root + url,
            "auto_toc_tree_section": "Contents",
        },
        True,
    )
    app.add_transform(AutoStructify)
