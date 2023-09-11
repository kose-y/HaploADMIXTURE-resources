using Random
using SnpArrays
using OpenADMIXTURE
using HaploADMIXTURE
using CUDA
using DelimitedFiles
function main()
    Ks = [2, 3, 4, 5, 6, 7, 8]
    snps = 100000
    nblocks = [snps ÷ 2]
    alphas = [0.1, 0.05, 0.02]
    seed = 95376
    for K in Ks
    for alpha in alphas
    for nblock in nblocks    
        println("running: K=$(K)")
	I = 1000
        J = nblock
        T = Float64
        T2 = Float32
        g = SnpLinAlg{T2}(SnpArray("g_5pops_1000indiv_$(snps)snps_$(seed)_$(alpha).bed"))#SnpArrays.datadir("EUR_subset.bed"))")))#"/home/kose/SKFR/1000genomes_new/SCOPE/misc/simulations/n1k_10k_5pops_final.bed"))#SnpArrays.datadir("EUR_subset.bed")));
        rng = MersenneTwister(3517)
        Q = 3
        d = HaploADMIXTURE.AdmixData2{T, T2}(I, J, K, Q, g; rng=rng)
        d_cu, g_cu = HaploADMIXTURE._cu_admixture_base(d, g, d.I, d.J)
#     d_cu, g_cu = nothing, nothing
        HaploADMIXTURE.init_em!(d, g, 5; d_cu=d_cu, g_cu=g_cu)
        HaploADMIXTURE.admixture_qn!(d, g, 300, 1e-7; d_cu=d_cu, g_cu=g_cu, mode=:ZAL, penalty=nothing)
        
        q = copy(d.q)
        p = copy(d.p)
        writedlm("cv/sim_5pops_1000indiv_$(snps)snps_$(seed)_$(alpha)_K$(K)_haplo_q.txt", q)
        writedlm("cv/sim_5pops_1000indiv_$(snps)snps_$(seed)_$(alpha)_K$(K)_haplo_p.txt", p)
    end
    end
    end
end
main()
