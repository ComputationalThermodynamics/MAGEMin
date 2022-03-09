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
sys.path.insert(0, os.path.abspath('.'))

from clang.cindex import Config

# 1. set the lib name, assuming it is in search path
Config.set_library_file('/usr/lib/llvm-6.0/lib/libclang.so.1')



# -- Project information -----------------------------------------------------

project = 'MAGEMin'
copyright = '2021, Riel, N., & Kaus, B.'
author = 'Riel, N., & Kaus, B.'

# The full version, including alpha/beta/rc tags
release = '2021'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
	'sphinx.ext.autodoc',
	'sphinx.ext.coverage',
	'sphinx.ext.napoleon',
	'sphinx_c_autodoc',
	'sphinx_c_autodoc.napoleon',
	'sphinx_c_autodoc.viewcode',
	'sphinx.ext.imgmath',
]

default_role = 'math'

imgmath_font_size   = 18  # for font size 14

#imgmath_dvipng_args = ['-gamma', '1.5', '-D', '300', '-bg', 'Transparent']
imgmath_image_format = 'svg'

# imgmath_latex = 'xelatex'
# imgmath_latex_args = ['--no-pdf']
# imgmath_latex_preamble = '\\usepackage{unicode-math}\\setmathfont{XITS Math}'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']




# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

html_theme_options = {
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# These folders are copied to the documentation's HTML output
html_static_path = ['_static']

html_css_files = [
    "css/custom.css"
]

# Theme options
html_theme_options = {
    # if we have a html_logo below, this shows /only/ the logo with no title text
    "logo_only": True,
    # Collapse navigation (False makes it tree-like)
    "collapse_navigation": False,
}

html_logo = "figs/logo_v3_text_white_bg.png"

c_autodoc_roots = ['../src/']

set_type_checking_flag = True
autodoc_member_order = 'bysource'
