"""WISE v3 driver (joint-conventions rebuild): build and inject only.

The standalone WISE survey is the frozen v2.1 script pipeline
(scripts/run_v2.sh); this driver exists for the joint v3 stage's WISE
tensor/injection inputs. ``freeze`` writes the operational v3 freeze
(joint split); the analysis stages are refused — the decision rule
lives in surveys/joint/joint.py.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from profile import PROFILE, write_freeze  # noqa: E402

from sglsurvey.cli import main  # noqa: E402

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "freeze":
        write_freeze()
    elif len(sys.argv) > 1 and sys.argv[1] in ("nulls", "completeness", "adjudicate", "report", "control"):
        raise SystemExit("wise v3 has no standalone decision rule — use surveys/joint/joint.py")
    else:
        main(PROFILE)
