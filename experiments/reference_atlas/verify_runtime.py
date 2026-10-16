"""Check a frozen model's runtime and export its public numeric parameters."""
import argparse
import hashlib
import json
from importlib.metadata import version
from pathlib import Path


def verify(model_path, output):
    import joblib
    import numpy as np
    from scipy.special import expit
    if hashlib.sha256(Path(model_path).read_bytes()).hexdigest() != "ae8673bd3a42a4220abcd35fed00f60c5e93c7ad69c9ae77cb2ae021ee3247d8":
        raise ValueError("Expected the original TTU model checkpoint")
    bundle = joblib.load(model_path)
    classifier = bundle["model"]
    scaler, logistic = classifier.steps[0][1], classifier.steps[1][1]
    if version("scikit-learn") != "1.7.2" or bundle["schema"] != "full-profile-transfer-v1":
        raise ValueError("Unexpected runtime/model schema")
    features = bundle["features"]
    if len(features) != 7 or logistic.classes_.tolist() != [0, 1]:
        raise ValueError("Unexpected features or classes")
    # Numeric implementation check only. Random feature vectors are not a
    # biological benchmark and are never used to tune any prediction rule.
    x = np.random.default_rng(20261016).normal(size=(64, 7))
    direct = expit(((x - scaler.mean_) / scaler.scale_) @ logistic.coef_[0] + logistic.intercept_[0])
    predicted = classifier.predict_proba(x)[:, 1]
    np.testing.assert_allclose(direct, predicted, rtol=1e-12, atol=1e-14)
    record = {"model_sha256": hashlib.sha256(Path(model_path).read_bytes()).hexdigest(),
        "software": {p: version(p) for p in ("numpy", "scipy", "scikit-learn", "biopython", "joblib", "pandas", "threadpoolctl")},
        "schema": bundle["schema"], "features": features, "threshold": bundle["threshold"],
        "scaler_mean": scaler.mean_.tolist(), "scaler_scale": scaler.scale_.tolist(),
        "logistic_coef": logistic.coef_[0].tolist(), "logistic_intercept": float(logistic.intercept_[0]),
        "numeric_equivalence_max_abs_error": float(np.max(np.abs(direct - predicted))),
        "numeric_equivalence_checked": True, "truth_used": False}
    Path(output).write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"numeric_equivalence_checked": True, "software": record["software"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    verify(args.model, args.output)
