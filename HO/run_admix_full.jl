using LinearAlgebra, Random, SnpArrays, LoopVectorization
using OpenADMIXTURE
using CUDA, CSV, DelimitedFiles
using StableRNGs

CUDA.allowscalar(true)
filename = "/mnt/Data2/HO/ho_final.bed"
for seed in [1234, 16962, 3517, 71000, 95376]
rng = StableRNG(seed)
d, _, _ = OpenADMIXTURE.run_admixture(filename, 12; T=Float32, use_gpu=true, rng=rng, admix_n_em_iters=3)
mkpath("adm")
writedlm( "adm/TGP_full.Q.12.$(seed).csv",  transpose(d.q), ',')
writedlm( "adm/TGP_full.P.12.$(seed).csv",  transpose(d.p), ',')
end
