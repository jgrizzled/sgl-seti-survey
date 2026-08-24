"""SPHEREx v2 driver: uv run python surveys/spherex/run.py <stage> [opts]."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from profile import PROFILE  # noqa: E402

from sglsurvey.v2.cli import main  # noqa: E402

if __name__ == "__main__":
    main(PROFILE)
