"""Checks for the reader-facing documentation expected in package code."""

import ast
from pathlib import Path


PACKAGE_DIR = Path(__file__).parents[1] / "xsuite_recycler"


def test_package_modules_and_top_level_callables_have_docstrings():
    """Require docs on modules, module functions, classes, and their methods."""
    missing = []
    for path in sorted(PACKAGE_DIR.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if ast.get_docstring(tree) is None:
            missing.append(f"{path.name}: module")

        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                if ast.get_docstring(node) is None:
                    missing.append(f"{path.name}:{node.lineno} {node.name}()")
            elif isinstance(node, ast.ClassDef):
                if ast.get_docstring(node) is None:
                    missing.append(f"{path.name}:{node.lineno} class {node.name}")
                for member in node.body:
                    if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        if ast.get_docstring(member) is None:
                            missing.append(
                                f"{path.name}:{member.lineno} "
                                f"{node.name}.{member.name}()"
                            )

    assert not missing, "Missing docstrings:\n" + "\n".join(missing)
