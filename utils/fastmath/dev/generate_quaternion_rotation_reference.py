# Rotation reference values for fastmath.quaternion (test/resources/quaternion/quaternion_rotation_reference.edn).
#
# Tools: Python scipy 1.18 (scipy.spatial.transform.Rotation) and pyquaternion 0.9.9 (/home/ts/penv, `uv run python`).
# scipy stores quaternions scalar-last (x, y, z, w), fastmath scalar-first (w, x, y, z) = (a, b, c, d): converted here.
# Domain: 60 random unit quaternions (numpy seed 12345) away from the gimbal lock and pairs for the interpolation.
# Conventions: intrinsic ZYX, angles [yaw pitch roll] in scipy, fastmath `[roll pitch yaw]` (to-euler, from-euler);
# intrinsic XYZ, [x y z] in both for to-angles and from-angles (R = Rx Ry Rz).
# Output: one map with vectors: :q (unit quaternions, scalar first, the sign as given by scipy `as_quat(canonical=True)`
# i.e. w >= 0), :euler [roll pitch yaw], :angles [x y z], :matrix (row-major 9 numbers), :axis-angle [angle ax ay az],
# :slerp-args [[q0 q1 t] ...] and :slerp (pyquaternion `Quaternion.slerp`, the shortest arc).
# Regeneration is byte-identical (fixed seed, repr of doubles).
import numpy as np
from scipy.spatial.transform import Rotation as R
from pyquaternion import Quaternion

rng = np.random.default_rng(12345)
rots = R.random(60, rng=rng)
q_xyzw = rots.as_quat(canonical=True)
qs = [[float(q[3]), float(q[0]), float(q[1]), float(q[2])] for q in q_xyzw]

zyx = rots.as_euler('ZYX')            # yaw, pitch, roll
euler = [[float(e[2]), float(e[1]), float(e[0])] for e in zyx]
xyz = rots.as_euler('XYZ')            # intrinsic x, y', z''
angles = [[float(e[0]), float(e[1]), float(e[2])] for e in xyz]
mats = [[float(x) for x in m.reshape(-1)] for m in rots.as_matrix()]
rv = rots.as_rotvec()
axis_angle = []
for v in rv:
    ang = float(np.linalg.norm(v))
    ax = v / ang
    axis_angle.append([ang, float(ax[0]), float(ax[1]), float(ax[2])])

slerp_args = []
slerp = []
for i in range(0, 58, 2):
    q0, q1 = qs[i], qs[i + 1]
    # a pair with a negative dot product is included on purpose: the sign of q1 is flipped for odd i
    if (i // 2) % 2 == 1:
        q1 = [-x for x in q1]
    for t in (0.1, 0.37, 0.5, 0.9):
        slerp_args.append([q0, q1, t])
        r = Quaternion.slerp(Quaternion(q0), Quaternion(q1), t)
        slerp.append([float(r[0]), float(r[1]), float(r[2]), float(r[3])])


def vec(v):
    return '[' + ' '.join(repr(float(x)) for x in v) + ']'


def tab(rows):
    return '[' + ' '.join(vec(r) for r in rows) + ']'


print('{')
print(':q ' + tab(qs))
print(':euler ' + tab(euler))
print(':angles ' + tab(angles))
print(':matrix ' + tab(mats))
print(':axis-angle ' + tab(axis_angle))
print(':slerp-args [' + ' '.join('[' + vec(a) + ' ' + vec(b) + ' ' + repr(t) + ']' for a, b, t in slerp_args) + ']')
print(':slerp ' + tab(slerp))
print('}')
