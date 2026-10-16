# Frozen full-profile match verifier

`model.joblib` is the original reference-trained checkpoint from the TTU
100-mixture experiment, retrieved without retraining for the UMD real-mixture
validation. SHA-256:

`ae8673bd3a42a4220abcd35fed00f60c5e93c7ad69c9ae77cb2ae021ee3247d8`

Training used Python 3.13.7, NumPy 2.3.3, scikit-learn 1.7.2, Biopython 1.85
and joblib 1.5.2. Features and prediction semantics are defined by
`experiments/reference_atlas/transfer.py`. The serialized bundle contains a
standardized logistic classifier, feature names, schema and reference-only
calibrated cutoff (0.94 up to floating-point representation).

Full [training/calibration provenance](../../results/transfer100_20261016/training.json)
and [original protocol](../../results/transfer100_20261016/PROTOCOL.json) are
preserved. Mixture truth was not used to fit this checkpoint. The learned score
is a candidate-match diagnostic, not a calibrated probability of biological
family presence. Apply it only with the documented patristic-distance feature
construction and nearest-family proposal; do not use it to rank distant families.

`numeric_parameters.json` preserves the scaler means/scales and logistic
coefficients as readable numbers. The UMD runtime check reconstructed scores
from these public parameters and matched the serialized classifier exactly on
64 deterministic numeric test vectors. This checks inference implementation,
not biological accuracy, and does not alter the model or its threshold.
