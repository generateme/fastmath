# Evaluation counts of scipy.optimize.toms748 (xtol 2e-12, rtol 8.9e-16) used by test/fastmath/solver/toms748_test.clj
# run: cd /home/ts/penv && uv run python <this file>
from scipy.optimize import toms748
import math
fs = {"cos x - x": (lambda x: math.cos(x) - x, 0, 1),
      "x^3 - 2": (lambda x: x**3 - 2, 0, 3),
      "tanh steep": (lambda x: math.tanh(100 * (x - .3)), 0, 1),
      "(x-1)^3": (lambda x: (x - 1)**3, 0, 3),
      "exp x - 1e10": (lambda x: math.exp(x) - 1e10, 0, 100),
      "atan x - 1": (lambda x: math.atan(x) - 1, 0, 10),
      "x^3 + x": (lambda x: x**3 + x, -1, 2),
      "1/x - 2": (lambda x: 1 / x - 2, .1, 5),
      "(x-1)^2 (x-3)": (lambda x: (x - 1)**2 * (x - 3), 2, 5),
      "x^5 - x - 1": (lambda x: x**5 - x - 1, 1, 2)}
for n, (f, a, b) in fs.items():
    r = toms748(f, a, b, xtol=2e-12, rtol=8.9e-16, full_output=True)
    print('["%s" %d]' % (n, r[1].function_calls))
