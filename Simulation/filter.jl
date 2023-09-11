using OpenADMIXTURE, SnpArrays, Random
using SKFR
using StableRNGs
using SnpArrays
prefix = ARGS[1]# SNPArrays.datadir("EUR_subset.bed")
X = SnpArray(prefix * ".bed")
I = size(X, 1)
n_pairs = parse(Int, ARGS[2]) ÷ 2
nclusters = parse(Int, ARGS[3]);

rng = StableRNG(7653)
IM = SKFR.get_imputed_matrix(X, nclusters; rng=rng, blocksize=2)

nblocks = [Int(x) for x in ([0.2, 0.15, 0.1, 0.05] .* n_pairs)]
_, selectedvecs = SKFR.sparsekmeans_path(IM, nblocks; iter=10);

for (nblock, selectedvec) in zip(nblocks, selectedvecs)
    aimlist_sorted = sort(selectedvec)
    des = "$(prefix)_$(nblock * 2)aims_ha"
    SnpArrays.filter(prefix, trues(I), aimlist_sorted; des=des)
end


IM = SKFR.get_imputed_matrix(X, nclusters; rng=rng, blocksize=1)

nblocks = [Int(x) for x in ([0.2, 0.15, 0.1, 0.05] .* n_pairs .* 2)]
_, selectedvecs = SKFR.sparsekmeans_path(IM, nblocks; iter=10);

for (nblock, selectedvec) in zip(nblocks, selectedvecs)
    aimlist_sorted = sort(selectedvec)
    des = "$(prefix)_$(nblock)aims_oa"
    SnpArrays.filter(prefix, trues(I), aimlist_sorted; des=des)
end
