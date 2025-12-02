import sys
import os
import shutil
import toml
import subprocess
import re
from pathlib import Path

FILE_PATH = os.path.dirname(__file__)

# Auto generate API documentation
def generate_api_docs():
    src_path = os.path.join(FILE_PATH, os.pardir, "TidalPy")
    out_path = os.path.join('API', 'generated')
    Path(out_path).mkdir(parents=True, exist_ok=True)

    subprocess.call([
        "sphinx-apidoc",
        "-o", str(out_path),
        str(src_path),
        "--force",
        "--implicit-namespaces",
        "--module-first",
        "--no-toc"
    ])
generate_api_docs()

# Basic configurations
sys.path.insert(0, os.path.abspath('../TidalPy'))
html_static_path = ["_static"]
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
html_logo = "_static/images/2025-11-28_Logo_2-4.svg"
pyproject_path = os.path.abspath(os.path.join(FILE_PATH, '..', 'pyproject.toml'))
with open(pyproject_path, 'r') as f:
    pyproject = toml.load(f)

project = pyproject['project']['name']
release = pyproject['project']['version']
author = 'Joe P. Renaud'

# Make a copy of the current change log and move it into docs so it can be included in the documentation.
src = os.path.abspath(os.path.join("..", "CHANGES.md"))
dst = os.path.abspath(os.path.join(FILE_PATH, "Changes.md"))
shutil.copyfile(src, dst)

src = os.path.abspath(os.path.join("..", "README.md"))
readme_file = os.path.abspath(os.path.join(FILE_PATH, "Overview", "Readme.md"))
shutil.copyfile(src, readme_file)

# For the readme we wan to decrement all headers by 1 so that they are loaded into the doctree correctly.
# Doing this here, rather than in the original readme.md file allows it to still serve as a stand alone document on 
# the repo's github and other locations.
# Read the copied file
with open(readme_file, 'r', encoding='utf-8') as f:
    content = f.read()

# Convert headers (remove one # from each header level)
# This pattern matches lines starting with 2+ hashes followed by a space
def reduce_header_level(match):
    hashes = match.group(1)
    rest = match.group(2)
    # Remove one hash (but keep at least one)
    new_hashes = hashes[1:] if len(hashes) > 1 else '#'
    return f"{new_hashes}{rest}"

# Apply the conversion
content = re.sub(r'^(#{2,})( .*)$', reduce_header_level, content, flags=re.MULTILINE)

# Write back to the file
with open(readme_file, 'w', encoding='utf-8') as f:
    f.write(content)

src = os.path.abspath(os.path.join("..", "LICENSE.md"))
dst = os.path.abspath(os.path.join(FILE_PATH, "Overview", "License.md"))
shutil.copyfile(src, dst)

src = os.path.abspath(os.path.join("..", "CONTRIBUTING.md"))
dst = os.path.abspath(os.path.join(FILE_PATH, "Overview", "Contributing.md"))
shutil.copyfile(src, dst)

src = os.path.abspath(os.path.join("..", "CODE_OF_CONDUCT.md"))
dst = os.path.abspath(os.path.join(FILE_PATH, "Overview", "CoC.md"))
shutil.copyfile(src, dst)

src = os.path.abspath(os.path.join("..", "NOTICE"))
dst = os.path.abspath(os.path.join(FILE_PATH, "Overview", "Notice.md"))
shutil.copyfile(src, dst)

# Copy all demos over
demo_src_dir = Path(os.path.join(FILE_PATH, os.pardir, "Demos"))
demo_dst_dir = Path(os.path.join(FILE_PATH, "Demos"))
demo_dst_dir.mkdir(exist_ok=True)
for notebook in demo_src_dir.glob("*.ipynb"):
    # New filename with underscores instead of spaces
    new_name = notebook.name.replace(" ", "_")
    dest_file = demo_dst_dir / new_name
    shutil.copy2(notebook, dest_file)

# Prep other theme and extension settings.
extensions = []

# Use Markdown files rather than rst (or in addition to)
extensions.append('myst_parser')
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
myst_enable_extensions = [
    "colon_fence",      # ::: fenced directives
    "deflist",          # definition lists
    "linkify",          # auto-detect URLs
    "smartquotes",      # nicer quotes
]

# Autodoc settings
extensions.append('sphinx.ext.autodoc')
extensions.append('sphinx.ext.viewcode')
extensions.append('sphinx.ext.autosummary')
autosummary_generate = True
autosummary_imported_members = True
autosummary_generate_recursive = True
autosummary_ignore_top = False
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "private-members": False,
    "show-inheritance": True,
}
napoleon_google_docstring = True
napoleon_numpy_docstring = True
# Support C++ autodocs
extensions.append('breathe')
breathe_default_project = "TidalPy"

# Jupyter notebook rendering
extensions.append('nbsphinx')
extensions.append('sphinx.ext.napoleon')
nbsphinx_allow_errors = True  # set True if you want docs to build even if notebooks fail
nbsphinx_execute = "auto"  # or "always"

# Copy code QOL button
extensions.append('sphinx_copybutton')
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True

# ReadTheDocs Themeing
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
}

# Looks like this has been deprecated. TODO: new analytics solution. 
# google_analytics_id = os.getenv("GOOGLE_ANALYTICS_ID", False)
# if google_analytics_id:
#     html_theme_options['analytics_id'] = google_analytics_id
#     html_theme_options['analytics_anonymize_ip'] = False



# Add custom CSS
def setup(app):
    app.add_css_file("custom.css")