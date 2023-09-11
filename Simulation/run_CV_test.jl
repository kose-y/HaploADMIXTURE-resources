using MLBase, StableRNGs
using SnpArrays, DelimitedFiles, HaploADMIXTURE, CUDA
function main()
I = 1000
nblocks = 50000
rng = StableRNG(100)
nfolds = 4
kf = collect(Kfold(rng, I, nfolds))
snps = 2 * nblocks
seed = 95376
alphas = [0.1, 0.05, 0.02]
for alpha in alphas
for K in 2:8
    q = readdlm("cv/sim_5pops_1000indiv_$(snps)snps_$(seed)_$(alpha)_K$(K)_haplo_q.txt")
    for fold in 1:nfolds
        p = readdlm("cv/train_5pops_1000indiv_$(snps)snps_$(seed)_$(alpha)_fold$(fold)_K$(K)_haplo_p.txt")
        println("running: K=$(K), fold=$(fold)")
        train = kf[fold]
        test = sort(setdiff(1:I, kf[fold]))
        I_test = length(test)
        J = nblocks
        T = Float64
        T2 = Float32
        g = SnpLinAlg{T2}(SnpArray("cv/test_5pops_1000indiv_$(snps)snps_$(seed)_$(alpha)_fold$(fold).bed"))#SnpArrays.datadir("EUR_subset.bed"))")))#"/home/kose/SKFR/1000genomes_new/SCOPE/misc/simulations/n1k_10k_5pops_final.bed"))#SnpArrays.datadir("EUR_subset.bed")));
        rng = StableRNG(3517)
        Q = 3
        d = HaploADMIXTURE.AdmixData2{T, T2}(I_test, J, K, Q, g; rng=rng)
        d_cu, g_cu = HaploADMIXTURE._cu_admixture_base(d, g, d.I, d.J)
#     d_cu, g_cu = nothing, nothing
        d.p .= p
        d.q .= @view(q[:, test])
        #HaploADMIXTURE.init_em!(d, g, 0; d_cu=d_cu, g_cu=g_cu, fix_p=true)
        HaploADMIXTURE.admixture_qn!(d, g, 300, 1e-7; d_cu=d_cu, g_cu=g_cu, mode=:ZAL, penalty=nothing, fix_p=true)
        q_test = copy(d.q)
        writedlm("cv/test_5pops_1000indiv_$(snps)snps_$(seed)_$(alpha)_fold$(fold)_K$(K)_haplo_q.txt", q_test)
    end
end
end
end
main()
