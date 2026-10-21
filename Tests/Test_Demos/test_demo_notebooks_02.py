"""Static checks of every committed demo and benchmark notebook, without running them: markdown free of control
characters (a tab where a LaTeX backslash-t was meant renders "\\times" as "imes"), relative links and their anchors
resolve, every code cell was executed and holds no error output, and no output shows a local path."""
import json
import re

import pytest

from demo_notebooks import NOTEBOOKS

LINK = re.compile(r"\]\(([^)\s]+)\)")
CONTROL = re.compile(r"[\x00-\x09\x0b-\x1f\x7f]")
# A user's home or a Windows drive path in an output: the notebook was run outside a clean checkout.
LOCAL_PATH = re.compile(r"[A-Za-z]:[\\/](Users|Documents and Settings)[\\/]|/(home|Users)/[A-Za-z]")


def _cells(path):
    return json.loads(path.read_text(encoding="utf-8"))["cells"]


def _slug(heading):
    """A heading's anchor as myst-parser and GitHub write it: lowercase, markup and punctuation dropped, spaces to -."""
    text = re.sub(r"[`*_]", "", heading.strip().lower())
    return re.sub(r"[^\w\- ]", "", text).replace(" ", "-")


def _anchors(path):
    if path.suffix == ".ipynb":
        lines = [line for cell in _cells(path) if cell["cell_type"] == "markdown"
                 for line in "".join(cell["source"]).splitlines()]
    else:
        lines = path.read_text(encoding="utf-8").splitlines()
    return {_slug(match.group(1)) for match in (re.match(r"#+\s+(.*)", line) for line in lines) if match}


@pytest.mark.parametrize("notebook_path", NOTEBOOKS, ids=[path.name for path in NOTEBOOKS])
def test_notebook_is_clean_without_running_it(notebook_path):
    problems = []
    for index, cell in enumerate(_cells(notebook_path)):
        source = "".join(cell["source"])
        if cell["cell_type"] == "markdown":
            if CONTROL.search(source):
                problems.append(f"cell {index}: control character in the markdown")
            for target in LINK.findall(source):
                if target.startswith(("http://", "https://", "mailto:", "attachment:")):
                    continue
                file_part, _, anchor = target.partition("#")
                linked = (notebook_path.parent / file_part).resolve() if file_part else notebook_path
                if not linked.exists():
                    problems.append(f"cell {index}: link to a missing file, {target}")
                elif anchor and linked.suffix in (".md", ".ipynb") and anchor not in _anchors(linked):
                    problems.append(f"cell {index}: link to a missing anchor, {target}")
        elif cell["cell_type"] == "code" and source.strip():
            if cell.get("execution_count") is None:
                problems.append(f"cell {index}: code cell not executed")
            for output in cell.get("outputs", []):
                if output.get("output_type") == "error":
                    problems.append(f"cell {index}: error output ({output.get('ename')})")
                text = "".join(output.get("text", "")) + "".join(output.get("data", {}).get("text/plain", ""))
                if LOCAL_PATH.search(text):
                    problems.append(f"cell {index}: local path in an output")
    assert not problems, "\n".join(problems)
