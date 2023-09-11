using Random
using SnpArrays
using OpenADMIXTURE
using HaploADMIXTURE
using CUDA
using DelimitedFiles
using StableRNGs
function main()
    Ks = 2:17
    nblocks = [385089 ÷ 2]
    seeds = [1234, 16962, 3517, 71000, 95376]
    for seed in seeds
    for K in Ks
    for nblock in nblocks    
        println("running: K=$(K)")
        println("seed: $(seed)")
	I = 1931
        J = nblock
        T = Float64
        T2 = Float32
        g = SnpLinAlg{T2}(SnpArray("/mnt/Data2/HO/ho_final_even.bed"))#SnpArrays.datadir("EUR_subset.bed"))")))#"/home/kose/SKFR/1000genomes_new/SCOPE/misc/simulations/n1k_10k_5pops_final.bed"))#SnpArrays.datadir("EUR_subset.bed")));
        rng = StableRNG(seed)
        Q = 3
        d = HaploADMIXTURE.AdmixData2{T, T2}(I, J, K, Q, g; rng=rng)
        d_cu, g_cu = HaploADMIXTURE._cu_admixture_base(d, g, d.I, d.J)
#     d_cu, g_cu = nothing, nothing
        HaploADMIXTURE.init_em!(d, g, 5; d_cu=d_cu, g_cu=g_cu)
        HaploADMIXTURE.admixture_qn!(d, g, 300, 1e-7; d_cu=d_cu, g_cu=g_cu, mode=:ZAL, penalty=nothing)
        
        q = copy(d.q)
        p = copy(d.p)
        writedlm("HO_1931_$(2 * nblock)_K$(K)_$(seed)_haplo_q.txt", q)
        writedlm("HO_1931_$(2 * nblock)_K$(K)_$(seed)_haplo_p.txt", p)
    end
    end
    end
end
main()
