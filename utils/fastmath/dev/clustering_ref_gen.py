"""Reference data for test/fastmath/ml/clustering_test.clj

Tools: Python 3.14, numpy, scikit-learn 1.9.0 (DBSCAN, KMeans, silhouette_samples,
NearestNeighbors) and kneed 0.9.0 (KneeLocator).
Run: cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/clustering_ref_gen.py
Output: test/resources/clustering/clustering_reference.edn (regeneration is byte-identical:
legacy RandomState streams, fixed order, repr() of doubles).

Why: fastmath.ml.clustering wraps Apache Commons Math (not tested here) and implements the
knee (Kneedle) radius selection, silhouette and wcss itself. Oracles:
- knee: kneed.KneeLocator(range(n), sorted k-distances, S, curve="convex",
  direction="increasing", online=False).knee_y, k-distance = distance to the k-th
  neighbour without the point itself (sklearn NearestNeighbors with n_neighbors = k + 1);
- dbscan: sklearn DBSCAN, min_samples = points + 1 (sklearn counts the point itself);
  labels, core sample indices (border points can be assigned differently by the two
  implementations, so the test compares core points and noise);
- silhouette: sklearn silhouette_samples / silhouette_score for the given labels;
- wcss: numpy, sum of squared distances to the mean of the cluster;
- kmeans: sklearn KMeans(n_init=10) labels on well separated data (inertia computed in numpy:
  KMeans.inertia_ is a threaded sum and differs in the last digit between runs).
"""
import os
import numpy as np
from sklearn.cluster import DBSCAN, KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.neighbors import NearestNeighbors
from kneed import KneeLocator

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "clustering", "clustering_reference.edn")


def blobs(seed, centers, sd, per):
    rs = np.random.RandomState(seed)
    return np.vstack([rs.normal(c, sd, size=(per, len(c))) for c in centers])


def datasets():
    d = {}
    d["blobs3-clean"] = blobs(1, [(0, 0), (6, 0), (3, 5)], 0.5, 30)
    rs = np.random.RandomState(11)
    d["blobs3"] = np.vstack([d["blobs3-clean"], rs.uniform([-4, -4], [10, 9], size=(6, 2))])
    d["blobs4-3d"] = blobs(3, [(0, 0, 0), (6, 0, 0), (0, 6, 0), (0, 0, 6)], 0.4, 25)
    rs = np.random.RandomState(2)
    a1, a2 = rs.uniform(0, 2 * np.pi, 60), rs.uniform(0, 2 * np.pi, 100)
    d["rings"] = np.vstack([np.c_[np.cos(a1), np.sin(a1)] + rs.normal(0, 0.08, (60, 2)),
                            3 * np.c_[np.cos(a2), np.sin(a2)] + rs.normal(0, 0.08, (100, 2))])
    rs = np.random.RandomState(5)
    t1, t2 = rs.uniform(0, np.pi, 50), rs.uniform(0, np.pi, 50)
    d["moons"] = np.vstack([np.c_[np.cos(t1), np.sin(t1)],
                            np.c_[1 - np.cos(t2), 0.5 - np.sin(t2)]]) + rs.normal(0, 0.05, (100, 2))
    d["uniform"] = np.random.RandomState(4).uniform(0, 1, (80, 2))
    d["expo"] = np.random.RandomState(6).exponential(1.0, (60, 1))
    return d


def k_distance_curve(x, k, metric):
    nn = NearestNeighbors(n_neighbors=k + 1, metric=metric).fit(x)
    dist, _ = nn.kneighbors(x)
    return np.sort(dist[:, k])


def knee_index(ys, s):
    kl = KneeLocator(list(range(len(ys))), ys, S=s, curve="convex",
                     direction="increasing", online=False)
    return None if kl.knee is None else int(kl.knee)


def knee(x, k, s, metric="euclidean"):
    """Returns [knee_y or None, robust]. A case is robust when the knee index does not
    change for S perturbed by 1e-9: k-distances of mutual nearest neighbours are exactly
    equal, which makes the Kneedle threshold test an analytic tie decided by the last
    digit of the distances (sklearn and Apache Commons Math differ there)."""
    ys = k_distance_curve(x, k, metric)
    i = knee_index(ys, s)
    lo, hi = (s * (1 - 1e-9), s * (1 + 1e-9)) if s > 0 else (-1e-9, 1e-9)
    robust = knee_index(ys, lo) == i and knee_index(ys, hi) == i
    return (None if i is None else float(ys[i])), robust


def edn(x):
    if x is None:
        return "nil"
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, (float, np.floating)):
        return repr(float(x))
    if isinstance(x, str):
        return '"' + x + '"'
    if isinstance(x, dict):
        return "{" + " ".join(":" + k + " " + edn(v) for k, v in x.items()) + "}"
    if isinstance(x, np.ndarray):
        return edn(x.tolist())
    return "[" + " ".join(edn(v) for v in x) + "]"


def wcss(x, labels, dist):
    total = 0.0
    for c in sorted(set(labels)):
        pts = x[labels == c]
        mu = pts.mean(axis=0)
        total += sum(dist(p, mu) ** 2 for p in pts)
    return float(total)


DIST = {"euclidean": lambda a, b: float(np.sqrt(((a - b) ** 2).sum())),
        "manhattan": lambda a, b: float(np.abs(a - b).sum()),
        "chebyshev": lambda a, b: float(np.abs(a - b).max())}


def main():
    ds = datasets()
    out = {"data": {k: v for k, v in ds.items()}}

    def knee_case(name, k, s, metric):
        eps, robust = knee(ds[name], k, s, metric)
        return {"data": name, "points": k, "sensitivity": s, "distance": metric,
                "eps": eps, "robust": robust}

    out["knee"] = [knee_case(name, k, s, "euclidean")
                   for name in ["blobs3", "blobs4-3d", "rings", "moons", "uniform", "expo"]
                   for k in [1, 2, 4, 8] for s in [0.0, 1.0, 3.0]]
    out["knee"] += [knee_case("blobs3", 4, 1.0, m) for m in ["manhattan", "chebyshev"]]

    dbs = [("blobs3", 0.8, 5, "euclidean"), ("blobs3", 1.0, 4, "manhattan"),
           ("blobs3", 0.5, 6, "chebyshev"), ("rings", 0.6, 4, "euclidean"),
           ("moons", 0.2, 5, "euclidean"), ("blobs4-3d", 1.0, 6, "chebyshev"),
           ("uniform", 0.08, 4, "euclidean")]
    out["dbscan"] = []
    for name, eps, ms, metric in dbs:
        m = DBSCAN(eps=eps, min_samples=ms, metric=metric).fit(ds[name])
        out["dbscan"].append({"data": name, "eps": eps, "points": ms - 1, "distance": metric,
                              "labels": m.labels_, "core": m.core_sample_indices_})

    rs = np.random.RandomState(7)
    cases = [("blobs3", KMeans(3, n_init=10, random_state=0).fit_predict(ds["blobs3"])),
             ("blobs3", KMeans(5, n_init=10, random_state=0).fit_predict(ds["blobs3"])),
             ("blobs4-3d", KMeans(4, n_init=10, random_state=0).fit_predict(ds["blobs4-3d"])),
             ("rings", rs.randint(0, 3, len(ds["rings"]))),
             ("uniform", np.r_[np.zeros(79, dtype=int), [1]])]
    out["silhouette"] = []
    for i, (name, labels) in enumerate(cases):
        for metric in (["euclidean", "manhattan", "chebyshev"] if i == 0 else ["euclidean"]):
            x = ds[name]
            out["silhouette"].append({
                "data": name, "labels": labels, "distance": metric,
                "mean": float(silhouette_score(x, labels, metric=metric)),
                "samples": silhouette_samples(x, labels, metric=metric),
                "wcss": wcss(x, labels, DIST[metric])})

    out["kmeans"] = []
    for name, k in [("blobs3-clean", 3), ("blobs4-3d", 4)]:
        km = KMeans(k, n_init=10, random_state=0).fit(ds[name])
        out["kmeans"].append({"data": name, "clusters": k, "labels": km.labels_,
                              "inertia": wcss(ds[name], km.labels_, DIST["euclidean"])})

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        f.write("{\n")
        keys = list(out.keys())
        for i, key in enumerate(keys):
            f.write(" :" + key + " " + edn(out[key]) + ("\n" if i == len(keys) - 1 else "\n"))
        f.write("}\n")


if __name__ == "__main__":
    main()
