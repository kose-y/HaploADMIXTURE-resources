using MLBase, StableRNGs
using SnpArrays, DelimitedFiles, HaploADMIXTURE, CUDA
function main()
I = 940
nblocks = 642950 ÷ 2
rng = StableRNG(100)
nfolds = 4
kf = collect(Kfold(rng, I, nfolds))
for K in 2:14
    p = readdlm("HGDP_940_642950_K$(K)_haplo_p.txt")
    q = readdlm("HGDP_940_642950_K$(K)_haplo_q.txt")
    for fold in 1:nfolds
        println("running: K=$(K), fold=$(fold)")
        train = kf[fold]
        test = setdiff(1:I, kf[fold])
        I_train = length(train)
        J = nblocks
        T = Float64
        T2 = Float32
        g = SnpLinAlg{T2}(SnpArray("train_fold$(fold).bed"))#SnpArrays.datadir("EUR_subset.bed"))")))#"/home/kose/SKFR/1000genomes_new/SCOPE/misc/simulations/n1k_10k_5pops_final.bed"))#SnpArrays.datadir("EUR_subset.bed")));
        rng = StableRNG(3517)
        Q = 3
        d = HaploADMIXTURE.AdmixData2{T, T2}(I_train, J, K, Q, g; rng=rng)
        d_cu, g_cu = HaploADMIXTURE._cu_admixture_base(d, g, d.I, d.J)
#     d_cu, g_cu = nothing, nothing
# HaploADMIXTURE.init_em!(d, g, 5; d_cu=d_cu, g_cu=g_cu)
        d.p .= p
        d.q .= @view(q[:, train])
        HaploADMIXTURE.admixture_qn!(d, g, 300, 1e-7; d_cu=d_cu, g_cu=g_cu, mode=:ZAL, penalty=nothing)
        q_train = copy(d.q)
        p_train = copy(d.p)
        writedlm("train_fold$(fold)_K$(K)_haplo_q.txt", q_train)
        writedlm("train_fold$(fold)_K$(K)_haplo_p.txt", p_train)
    end
end

end
main()
