using Random
using StableRNGs
using SnpArrays
using OpenADMIXTURE
using HaploADMIXTURE
using CUDA
using DelimitedFiles
function main()
    nblocks = [50000, 40000, 30000, 20000, 10000, 5000, 2500]
    for nblock in nblocks    
        I = 1718
        J = nblock
        K = 8
        T = Float64
        T2 = Float64
        g = SnpLinAlg{T2}(SnpArray(SnpArrays.datadir("/home/kose/HaploADMIXTURE/TGP/TGP_8_$(2 * nblock)aims.bed")))#SnpArrays.datadir("EUR_subset.bed"))")))#"/home/kose/SKFR/1000genomes_new/SCOPE/misc/simulations/n1k_10k_5pops_final.bed"))#SnpArrays.datadir("EUR_subset.bed")));
        rng = StableRNG(3517)
        Q = 3
        d = HaploADMIXTURE.AdmixData2{T, T2}(I, J, K, Q, g; rng=Random.GLOBAL_RNG)
        d_cu, g_cu = HaploADMIXTURE._cu_admixture_base(d, g, d.I, d.J)
#     d_cu, g_cu = nothing, nothing
        HaploADMIXTURE.init_em!(d, g, 5; d_cu=d_cu, g_cu=g_cu)
        HaploADMIXTURE.admixture_qn!(d, g, 200, 1e-7; d_cu=d_cu, g_cu=g_cu, mode=:ZAL, penalty=true)

        writedlm("TGP_1718_$(2 * nblock)_haplo_q_pen.txt", d.q)
        writedlm("TGP_1718_$(2 * nblock)_haplo_p_pen.txt", d.p)
    end
end
main()
