"""``python -m atmospheric_data <command> config/<case>.yaml`` entry point."""
import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
