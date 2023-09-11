using MLBase, StableRNGs
using SnpArrays, DelimitedFiles, HaploADMIXTURE, CUDA
function main()
I = 1931
nblocks = 385088 ÷ 2
rng = StableRNG(100)
nfolds = 4
kf = collect(Kfold(rng, I, nfolds))
for K in 2:17
    q = readdlm("HO_1931_385088_K$(K)_haplo_q.txt")
    for fold in 1:nfolds
        p = readdlm("train_fold$(fold)_K$(K)_haplo_p.txt")
        println("running: K=$(K), fold=$(fold)")
        train = kf[fold]
        test = sort(setdiff(1:I, kf[fold]))
        I_test = length(test)
        J = nblocks
        T = Float64
        T2 = Float32
        g = SnpLinAlg{T2}(SnpArray("test_fold$(fold).bed"))#SnpArrays.datadir("EUR_subset.bed"))")))#"/home/kose/SKFR/1000genomes_new/SCOPE/misc/simulations/n1k_10k_5pops_final.bed"))#SnpArrays.datadir("EUR_subset.bed")));
        rng = StableRNG(3517)
        Q = 3
        d = HaploADMIXTURE.AdmixData2{T, T2}(I_test, J, K, Q, g; rng=rng)
        d_cu, g_cu = HaploADMIXTURE._cu_admixture_base(d, g, d.I, d.J)
#     d_cu, g_cu = nothing, nothing
        #HaploADMIXTURE.init_em!(d, g, 5; d_cu=d_cu, g_cu=g_cu)
        d.p .= p
        d.q .= @view(q[:, test])
        HaploADMIXTURE.admixture_qn!(d, g, 300, 1e-7; d_cu=d_cu, g_cu=g_cu, mode=:ZAL, penalty=nothing, fix_p=true)
        q_test = copy(d.q)
        writedlm("test_fold$(fold)_K$(K)_haplo_q.txt", q_test)
    end
end

end
main()
