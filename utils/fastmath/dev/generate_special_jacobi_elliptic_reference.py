# Reference values for the Jacobi and elliptic fixes of fastmath.special (upstream group 5):
#
# - :am          Jacobi amplitude outside 0 <= m <= 1: the root of ellipf(phi, m) = u for m < 0; atan2(sn, cn) for m > 1
#                (bounded, the mpmath convention);
# - :inverse     the 12 inverse Jacobi functions at x = double(f(u)), the reference u is the root of f(u) = x on the
#                principal interval of the function (DLMF 22.15): [-K, K] for sn, sc, sd; [0, 2K] where f takes both signs
#                there; [0, K] for an odd or a positive-only f (negative x of an odd f from -u). K is the real quarter
#                period, K(1/m)/sqrt(m) for m > 1;
# - :complete    K, E, D = (K - E)/m and PI(0.5|m) at large negative m;
# - :near-one    sn, cn, dn, am next to m = 1 (am = j pi + atan(sn/cn) on the j-th half period, gd(u) at m = 1);
# - :negative-m  sn, cn, dn at large negative m (ill-conditioned: the phase is u sqrt(-m); 400 digits).
#
# Tool: Python mpmath 1.3.0, 50 decimal digits, at the exact double inputs. Run in the /home/ts/penv environment:
#
#   cd /home/ts/clojure/fastmath && /home/ts/penv/.venv/bin/python utils/fastmath/dev/generate_special_jacobi_elliptic_reference.py
#
# Output: test/resources/special/jacobi_elliptic_reference.edn

import os
import mpmath as mp

mp.mp.dps = 50


def fmt(x):
    x = float(x)
    if x == float('inf'):
        return '##Inf'
    if x == float('-inf'):
        return '##-Inf'
    return repr(x)


def vec(xs):
    return '[' + ' '.join(fmt(x) for x in xs) + ']'


def strvec(xs):
    return '[' + ' '.join('"%s"' % x for x in xs) + ']'


def jf(name, u, m):
    return mp.re(mp.ellipfun(name, u, m=m))


f = jf


def quarter(m):
    return mp.ellipk(1 / m) / mp.sqrt(m) if m > 1 else mp.ellipk(m)


def am_root(u, m):
    u, m = mp.mpf(u), mp.mpf(m)
    return mp.findroot(lambda p: mp.ellipf(p, m) - u, mp.pi * u / (2 * mp.ellipk(m)))


def am_half_period(u, m):
    # am - j pi in [-pi/2, pi/2] on [(2j - 1)K, (2j + 1)K]: am = j pi + atan(sn/cn) (Newton on ellipf fails near m = 1)
    u, m = mp.mpf(u), mp.mpf(m)
    K = mp.ellipk(m)
    j = mp.floor((u + K) / (2 * K))
    return j * mp.pi + mp.atan(jf('sn', u, m) / jf('cn', u, m))


# :am
au, am_, ar = [], [], []
for m in [-100.0, -2.0, -0.5, -1e-3]:
    for u in [-3.0, 0.7, 3.0, 10.0]:
        au.append(u); am_.append(m); ar.append(am_root(u, m))
for m in [1.5, 3.0, 10.0]:
    for u in [-3.0, 0.7, 3.0, 10.0]:
        au.append(u); am_.append(m)
        ar.append(mp.atan2(jf('sn', mp.mpf(u), mp.mpf(m)), jf('cn', mp.mpf(u), mp.mpf(m))))

# :inverse, principal interval class per function: 'odd' [-K, K]; 'both' [0, 2K]; 'half' [0, K] (+ -u for odd)
CLASS_LOW = {'sn': 'odd', 'sc': 'odd', 'sd': 'odd', 'cn': 'both', 'cs': 'both', 'cd': 'both', 'dc': 'both',
             'nc': 'both', 'ds': 'half-odd', 'ns': 'half-odd', 'dn': 'half', 'nd': 'half'}
CLASS_HIGH = {'sn': 'odd', 'sc': 'odd', 'sd': 'odd', 'dn': 'both', 'cd': 'both', 'ds': 'both', 'dc': 'both',
              'nd': 'both', 'cs': 'half-odd', 'ns': 'half-odd', 'cn': 'half', 'nc': 'half'}
FRACTIONS = [0.05, 0.2, 0.45, 0.7, 0.95]

inf_f, inf_x, inf_m, inf_r = [], [], [], []
for m in [-2.0, -0.5, 0.3, 0.8, 1.5, 3.0]:
    mm = mp.mpf(m)
    K = quarter(mm)
    for name in sorted(CLASS_LOW):
        cls = (CLASS_HIGH if m > 1 else CLASS_LOW)[name]
        if cls == 'odd':
            us = [(2 * t - 1) * K for t in FRACTIONS]
        elif cls == 'both':
            us = [2 * t * K for t in FRACTIONS]
        elif cls == 'half-odd':
            us = [t * K for t in FRACTIONS] + [-t * K for t in FRACTIONS]
        else:
            us = [t * K for t in FRACTIONS]
        for u in us:
            x = mp.mpf(float(f(name, u, mm)))
            r = mp.findroot(lambda v: f(name, v, mm) - x, u)
            assert abs(f(name, r, mm) - x) <= mp.mpf(1e-14) * max(1, abs(x)), (name, m, u)
            assert abs(r - u) < K / 10, (name, m, u, r)
            inf_f.append(name); inf_x.append(x); inf_m.append(m); inf_r.append(r)

# :complete
cm, cK, cE, cD, cP = [], [], [], [], []
for m in [-9.0, -10.0, -1e3, -1e6, -1e10, -1e12, -1e15, -1e16, -1e20, -1e300]:
    mm = mp.mpf(m)
    k, e = mp.ellipk(mm), mp.ellipe(mm)
    cm.append(m); cK.append(k); cE.append(e); cD.append((k - e) / mm); cP.append(mp.ellippi(mp.mpf(0.5), mm))

# :near-one
nu, nm, nsn, ncn, ndn, nam = [], [], [], [], [], []
for m in [1 - 1e-8, 1 - 1e-10, 1 - 1e-12, 1 - 1e-15, 1.0]:
    mm = mp.mpf(m)
    for u in [0.5, 5.0, 20.0, 30.0, 50.0]:
        uu = mp.mpf(u)
        nu.append(u); nm.append(m)
        nsn.append(jf('sn', uu, mm)); ncn.append(jf('cn', uu, mm)); ndn.append(jf('dn', uu, mm))
        nam.append(2 * mp.atan(mp.tanh(uu / 2)) if m == 1.0 else am_half_period(u, m))

# :negative-m
gu, gm, gsn, gcn, gdn = [], [], [], [], []
for m in [-1e3, -1e10, -1e17, -1e300]:
    mm = mp.mpf(m)
    for u in [1e-3, 0.3, 3.0]:
        uu = mp.mpf(u)
        gu.append(u); gm.append(m)
        with mp.workdps(400):  # the phase u sqrt(-m) reaches 3e150: 50 digits give a wrong dn < 0 at m = -1e300
            gsn.append(+jf('sn', uu, mm)); gcn.append(+jf('cn', uu, mm)); gdn.append(+jf('dn', uu, mm))

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                   'test', 'resources', 'special', 'jacobi_elliptic_reference.edn')
with open(out, 'w') as fo:
    fo.write('{:am {:u %s :m %s :ref %s}\n' % (vec(au), vec(am_), vec(ar)))
    fo.write(' :inverse {:family %s :x %s :m %s :ref %s}\n' % (strvec(inf_f), vec(inf_x), vec(inf_m), vec(inf_r)))
    fo.write(' :complete {:m %s :K %s :E %s :D %s :PI %s}\n' % (vec(cm), vec(cK), vec(cE), vec(cD), vec(cP)))
    fo.write(' :near-one {:u %s :m %s :sn %s :cn %s :dn %s :am %s}\n' % (vec(nu), vec(nm), vec(nsn), vec(ncn),
                                                                          vec(ndn), vec(nam)))
    fo.write(' :negative-m {:u %s :m %s :sn %s :cn %s :dn %s}}\n' % (vec(gu), vec(gm), vec(gsn), vec(gcn), vec(gdn)))
print('am', len(au), 'inverse', len(inf_f), 'complete', len(cm), 'near-one', len(nu), 'negative-m', len(gu))
