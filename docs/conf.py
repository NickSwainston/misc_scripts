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
import logging

logger = logging.getLogger(__name__)

DOC_SOURCES_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT_DIR = os.path.abspath(os.path.join(DOC_SOURCES_DIR, ".."))

# insert
sys.path.insert(0, DOC_SOURCES_DIR)
sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'misc_scripts'
copyright = '2022, Nick Swainston'
author = 'Nick Swainston'

# The full version, including alpha/beta/rc tags
release = '1.0'

# Fix from https://github.com/readthedocs/readthedocs.org/issues/1846
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    # from convert_readme import convert_md
    # convert_md()

    render_examples = False

    # hack for lacking git-lfs support on rtd
    import git_lfs
    from urllib.error import HTTPError

    _fetch_urls = git_lfs.fetch_urls

    def _patched_fetch_urls(lfs_url, oid_list):
        """Hack git_lfs library that sometimes makes too big requests"""
        objects = []

        try:
            objects.extend(_fetch_urls(lfs_url, oid_list))
        except HTTPError as err:
            if err.code != 413:
                raise
            logger.error("LFS: request entity too large, splitting in half")
            objects.extend(_patched_fetch_urls(lfs_url, oid_list[:len(oid_list) // 2]))
            objects.extend(_patched_fetch_urls(lfs_url, oid_list[len(oid_list) // 2:]))

        return objects

    git_lfs.fetch_urls = _patched_fetch_urls
    git_lfs.fetch(PROJECT_ROOT_DIR)

else:
    render_examples = True


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'numpydoc',
    'myst_parser'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']