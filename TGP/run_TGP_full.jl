using Random
using SnpArrays
using OpenADMIXTURE
using HaploADMIXTURE
using CUDA
using DelimitedFiles
function main()
    nblocks = [1854622 ÷ 2]
    for nblock in nblocks    
        I = 1718
        J = nblock
        K = 8
        T = Float64
        T2 = Float32
        g = SnpLinAlg{T2}(SnpArray("/mnt/Data2/TGP/TGP_1718.bed"))#SnpArrays.datadir("EUR_subset.bed"))")))#"/home/kose/SKFR/1000genomes_new/SCOPE/misc/simulations/n1k_10k_5pops_final.bed"))#SnpArrays.datadir("EUR_subset.bed")));
        rng = MersenneTwister(3517)
        Q = 3
        d = HaploADMIXTURE.AdmixData2{T, T2}(I, J, K, Q, g; rng=rng)
        d_cu, g_cu = HaploADMIXTURE._cu_admixture_base(d, g, d.I, d.J)
        HaploADMIXTURE.init_em!(d, g, 5; d_cu=d_cu, g_cu=g_cu)
        HaploADMIXTURE.admixture_qn!(d, g, 300, 1e-7; d_cu=d_cu, g_cu=g_cu, mode=:ZAL)

        writedlm("TGP_1718_$(2 * nblock)_haplo_q_32.txt", d.q)
        writedlm("TGP_1718_$(2 * nblock)_haplo_p_32.txt", d.p)
    end
end
main()
