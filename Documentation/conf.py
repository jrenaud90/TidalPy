import sys
import os
import shutil
sys.path.insert(0, os.path.abspath('../TidalPy'))
import toml

html_static_path = ["_static"]
pyproject_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'pyproject.toml'))
with open(pyproject_path, 'r') as f:
    pyproject = toml.load(f)

project = pyproject['project']['name']
release = pyproject['project']['version']
author = 'Joe P. Renaud'

# Make a copy of the current change log and move it into docs so it can be included in the documentation.
src = os.path.abspath(os.path.join("..", "CHANGES.md"))
dst = os.path.abspath(os.path.join(os.path.dirname(__file__), "Changes.md"))
shutil.copyfile(src, dst)

src = os.path.abspath(os.path.join("..", "README.md"))
dst = os.path.abspath(os.path.join(os.path.dirname(__file__), "Readme.md"))
shutil.copyfile(src, dst)

src = os.path.abspath(os.path.join("..", "LICENSE.md"))
dst = os.path.abspath(os.path.join(os.path.dirname(__file__), "License.md"))
shutil.copyfile(src, dst)

src = os.path.abspath(os.path.join("..", "NOTICE"))
dst = os.path.abspath(os.path.join(os.path.dirname(__file__), "Notice.md"))
shutil.copyfile(src, dst)

extensions = [
    'myst_parser',
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_logo = "_static/images/2025-11-28_Logo_2-4.svg"
