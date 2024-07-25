# tvCoef 0.2.1

* deprecated functions removed.

* correction in `bp_lm()` of the use of the intercept.

* correction of the use of `bp_lm()` in `rmse_prev`.

# tvCoef 0.2.0

* function `get_coeff_plot()` removed to avoid using `ggplot2` and `patchwork`.

* doc improvements.

* `prevision` attribute (for example in `oos_prev()`) renamed to `forecast`.

* renaming of exported S3 classes : `"bp.lms"`, `"piecereg"` and `"hansen.test"` to `"bp_lms"`, `"piecereg"` and `"hansen_test"` respectively.
