"""Label-free end-to-end PyEuk check, run INSIDE the container.

Mirrors /home/anton/pyeuk-bench/lofreq_arm/head_unsup.py exactly, with
find_clusters(matrix, None) -- no gold file passed to the clusterer.
"""
import os, sys, itertools, importlib.metadata
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine
from cyclospora_pyeuk.clustering import CyclosporaClusterFinder
import cyclospora_pyeuk

print("python      :", sys.version.split()[0])
print("pyeuk       :", cyclospora_pyeuk.__version__, "at", os.path.dirname(cyclospora_pyeuk.__file__))
for p in ("numpy", "scipy", "pandas", "scikit-learn", "numba", "pysam"):
    print("%-12s: %s" % (p, importlib.metadata.version(p)))

sheet, gold_path, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
gold = pd.read_csv(gold_path, sep="\t")
gmap = dict(zip(gold.Seq_ID, gold.Cluster_alias))
df = pd.read_csv(sheet, sep="\t")

mat = PyEukDistanceEngine(epsilon=0.3072).compute_revised_wibs_matrix(df)
CyclosporaClusterFinder(stringency=95.0, robust=True).find_clusters(mat, None, output_dir=out)

cl = None
for fn in sorted(os.listdir(out)):
    if "RESULTING_CLUSTERS" in fn:
        cl = pd.read_csv(os.path.join(out, fn), sep="\t")
assert cl is not None, "no RESULTING_CLUSTERS file written"
cl = cl[cl.Seq_ID.isin(gmap)]
yt = [gmap[i] for i in cl.Seq_ID]
yp = list(cl.Assigned_cluster)
n = len(yt)
TP = FP = FN = TN = 0
for i, j in itertools.combinations(range(n), 2):
    st, sp = yt[i] == yt[j], yp[i] == yp[j]
    TP += st and sp; FP += (not st) and sp
    FN += st and not sp; TN += (not st) and (not sp)
ari = adjusted_rand_score(yt, yp)
k = len(set(yp))
print("LABEL-FREE n=%d ARI=%.4f sens=%.4f spec=%.4f k=%d sizes=%s"
      % (n, ari, TP / (TP + FN), TN / (TN + FP), k, pd.Series(yp).value_counts().to_dict()))
ok = (k == 2 and ari >= 0.94)
print("VERDICT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
