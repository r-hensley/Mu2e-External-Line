"""Run the Recycler filamentation example with ``python -m``."""

try:
    from .run import main
except ImportError:  # Support direct execution with ``python __main__.py``.
    from run import main


if __name__ == "__main__":
    raise SystemExit(main())
