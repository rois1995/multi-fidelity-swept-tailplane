"""Reproduce the Table 2 fit interval without altering the experimental data."""
import numpy as np
from scipy.stats import t
from load import WT

w = WT()
selected = (w['Delta'] == 0) & (np.abs(w['AoA']) <= 10)
alpha = w['AoA'][selected]
# Match the unweighted slope-ratio estimate used for the published coordinate.
slopes, half_widths = [], []
for coefficient in ('CMy', 'CFz'):
    fit, covariance = np.polyfit(alpha, w[coefficient][selected], 1, cov=True)
    slopes.append(fit[0])
    half_widths.append(t.ppf(0.975, len(alpha) - 2) * np.sqrt(covariance[0, 0]))
moment_slope, force_slope = slopes
ratio = moment_slope / force_slope
chord = 0.679183673469
coordinate = 0.8008 - chord * ratio
# Propagate the 95% slope-fit half-widths assuming independent slope estimates.
interval = chord * abs(ratio) * np.linalg.norm(np.array(half_widths) / slopes)
assert np.isfinite(interval) and interval > 0
assert f'{coordinate:.4f}' == '0.6426'
assert f'{interval:.4f}' == '0.0011'
print(f'N={len(alpha)}; delta=0 deg; |alpha|<=10 deg; unweighted linear fits')
print(f'x_AC={coordinate:.10f} m; propagated fit interval=+/-{interval:.10f} m')
print('95% slope-fit intervals propagated without slope covariance; not a total measurement bound.')
