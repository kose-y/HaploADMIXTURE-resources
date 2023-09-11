using Random
using SnpArrays
using OpenADMIXTURE
using HaploADMIXTURE
using CUDA
using DelimitedFiles
function main()
    nblocks = [642951 ÷ 2]
    for nblock in nblocks    
        I = 940
        J = nblock
        K = 10
        T = Float64
        T2 = Float64
        g = SnpLinAlg{T2}(SnpArray("/mnt/Data2/HGDP/HGDP_940.bed"))#SnpArrays.datadir("EUR_subset.bed"))")))#"/home/kose/SKFR/1000genomes_new/SCOPE/misc/simulations/n1k_10k_5pops_final.bed"))#SnpArrays.datadir("EUR_subset.bed")));
        rng = MersenneTwister(3517)
        Q = 3
        d = HaploADMIXTURE.AdmixData2{T, T2}(I, J, K, Q, g; rng=rng)
        d_cu, g_cu = HaploADMIXTURE._cu_admixture_base(d, g, d.I, d.J)
#     d_cu, g_cu = nothing, nothing
        HaploADMIXTURE.init_em!(d, g, 5; d_cu=d_cu, g_cu=g_cu)
        HaploADMIXTURE.admixture_qn!(d, g, 300, 1e-7; d_cu=d_cu, g_cu=g_cu, mode=:ZAL)

        writedlm("HGDP_940_$(2 * nblock)_haplo_q.txt", d.q)
        writedlm("HGDP_940_$(2 * nblock)_haplo_p.txt", d.p)
    end
end
main()
