using OpenADMIXTURE, SnpArrays, Random
using SKFR, StableRNGs

using SnpArrays
filename = "/mnt/Data2/TGP/TGP_1718.bed"# SNPArrays.datadir("EUR_subset.bed")
X = SnpArray(filename)
I = size(X, 1)
nclusters = 7;

rng = MersenneTwister(7653)
IM = SKFR.get_imputed_matrix(X, nclusters; rng=rng, blocksize=2)

nblocks = [100000, 80000, 60000, 40000, 20000, 10000, 5000] .÷ 2
_, selectedvecs = SKFR.sparsekmeans_path(IM, nblocks; iter=10);

for (nblock, selectedvec) in zip(nblocks, selectedvecs)
    aimlist_sorted = sort(selectedvec)
    des = "TGP_$(nclusters)_$(nblock * 2)aims"
    SnpArrays.filter(filename[1:end-4], trues(I), aimlist_sorted; des=des)
end


