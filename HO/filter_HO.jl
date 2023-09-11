using OpenADMIXTURE, SnpArrays, Random
using SKFR

using SnpArrays
using StableRNGs
filename = "/mnt/Data2/HO/ho_final.bed"# SNPArrays.datadir("EUR_subset.bed")
X = SnpArray(filename)
I = size(X, 1)
nclusters = 12;

rng = StableRNG(7653)
IM = SKFR.get_imputed_matrix(X, nclusters; rng=rng, blocksize=1)

nblocks = [50000, 40000, 30000, 20000, 10000, 5000, 2500] .* 2
_, selectedvecs = SKFR.sparsekmeans_path(IM, nblocks; iter=10);

for (nblock, selectedvec) in zip(nblocks, selectedvecs)
    aimlist_sorted = sort(selectedvec)
    des = "adm/HO_$(nclusters)_$(nblock)aims"
    SnpArrays.filter(filename[1:end-4], trues(I), aimlist_sorted; des=des)
end


