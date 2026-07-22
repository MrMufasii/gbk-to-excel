#!/usr/bin/env python3
"""Execute the Aas class-specificity experiment without the optional IQ-TREE stage."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import aas_module_strengthening as analysis

analysis.build_tree = lambda *args, **kwargs: {
    "status": "skipped_in_fast_profile_run",
    "nearest_classes": {},
}
analysis.run_domain_specificity(Path("results/aas_profile_fast"))
