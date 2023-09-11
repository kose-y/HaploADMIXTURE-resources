using StatsBase
using Statistics
using SnpArrays
using NaNMath
using Random
using Missings
using LinearAlgebra
using Clustering
using DelimitedFiles
using SKFR
using StableRNGs
function shuffle(srcfilenm; desfilenm = split(srcfilenm, ".bed")[1] * ".shuffled.bed", rng=Random.GLOBAL_RNG)
    s_in = SnpArray(srcfilenm)
    m, n = size(s_in)
    s_out = SnpArray(desfilenm, m, n)
    
    @inbounds for j in 1:n
        order = randperm(rng, m)
        s_out[:, j] .= s_in[order, j]
    end
end




function run_SKFR(data_path, clusters, sparsities)

    # Run SKFR
    classes = clusters
    rng = StableRNG(1696)
    for i in 1:3
        shuffle(data_path; 
            desfilenm=split(data_path, ".bed")[1] * "_shuffled_$(i).bed", rng=rng)
    end

    rng = StableRNG(1696)
    classes = 1:10
    WSSs = zeros(Float64, length(classes))
    for sparsity in sparsities
    X = SnpArray(data_path)
    n, p = size(X)
    println("Sparsity: $(sparsity)")
    ISM_prev = nothing
    for i in reverse(classes)
        ISM = SKFR.ImputedSnpMatrix{Float64}(X, i; rng=rng, blocksize=2)
        if ISM_prev !== nothing 
            to_exclude = argmin(ISM_prev.members)
            if to_exclude != 1
                ISM.centers_tmp[:, 1:to_exclude-1] .= ISM_prev.centers_tmp[:, 1:to_exclude-1]
                ISM.centers[:, 1:to_exclude-1] .= ISM_prev.centers[:, 1:to_exclude-1]
            end
            if to_exclude != size(ISM_prev.centers_tmp, 2)
                ISM.centers_tmp[:, to_exclude:end] .= ISM_prev.centers_tmp[:, to_exclude+1:end]
                ISM.centers[:, to_exclude:end] .= ISM_prev.centers[:, to_exclude+1:end]
            end
            SKFR.get_distances_to_center!(ISM)
            SKFR.get_clusters!(ISM)
        end
        _, _, _, _WSS, _TSS = SKFR.sparsekmeans1(ISM, sparsity ÷ 2)
        WSSs[i] = sum(_WSS)
        ISM_prev = ISM
    end

    WSSs_shuffled = zeros(Float64, length(classes), 3)
    for rep in 1:3
        data_path_s = split(data_path, ".bed")[1] * "_shuffled_$(rep).bed"
        X = SnpArray(data_path_s, n)
            #Xtrue = convert(Matrix{Float64}, X, model=ADDITIVE_MODEL, center=false, scale=false)
        rng = StableRNG(1696)
        ISM_prev = nothing
        for i in reverse(classes)
            ISM = SKFR.ImputedSnpMatrix{Float64}(X, i; rng=rng, blocksize=2)
            if ISM_prev !== nothing 
                to_exclude = argmin(ISM_prev.members)
                if to_exclude != 1
                    ISM.centers_tmp[:, 1:to_exclude-1] .= ISM_prev.centers_tmp[:, 1:to_exclude-1]
                    ISM.centers[:, 1:to_exclude-1] .= ISM_prev.centers[:, 1:to_exclude-1]
                end
                if to_exclude != size(ISM_prev.centers_tmp, 2)
                    ISM.centers_tmp[:, to_exclude:end] .= ISM_prev.centers_tmp[:, to_exclude+1:end]
                    ISM.centers[:, to_exclude:end] .= ISM_prev.centers[:, to_exclude+1:end]
                end
                SKFR.get_distances_to_center!(ISM)
                SKFR.get_clusters!(ISM)
            end
            #Random.seed!(1696)
            _, _, _, _WSS, _TSS = SKFR.sparsekmeans1(ISM, sparsity ÷ 2)
            WSSs_shuffled[i, rep] =  sum(_WSS)
            ISM_prev = ISM
        end
    end
    r = Float64[]
    for i in 1:length(classes)
        append!(r, mean(log.(WSSs_shuffled[i,:])) - log.(WSSs[i]))
    end
    sd = std(log.(WSSs_shuffled), dims=2)
    r_minus_sd = r .- sd * sqrt(1 + 1/3)
    for i in 1:length(classes)-1
        if r[i] > r_minus_sd[i+1]
            println("Selected K : $i")
            break
        end
        if i == length(classes)-1
            println("Selected K: FULL($(i+1))")
        end
    end
    println(WSSs)
    println(WSSs_shuffled)
    println(r)
    println(r_minus_sd)
    println("END $(sparsity). ")
    flush(stdout)
    end

end


function run_case(n, p, k, prefixn, prefixp)
    data_path = "g_$(k)pops_$(n)indiv_$(p)snps_95376_0.1.bed"
    #data_path = "n$(prefixn)_$(prefixp)_$(k)pops_gamma001.bed"
    clusters = k

    sparsities = convert(Vector{Int}, [1.00, 0.50, 0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001] .* p)
    

    run_SKFR(data_path, clusters, sparsities);
    
    #f = open("rslts/clusters_sim_n$(prefixn)_$(prefixp)_$(k)pops_path", "w")
   # writedlm(f, clusters)
    #close(f)
#
    #for (v, s) in zip(snps, sparsities)
    #    f = open("rslts/aims_sim_n$(prefixn)_$(prefixp)_$(k)pops_path_$(s)", "w")
    #    writedlm(f, v)
    #    close(f)
    #    aims = sort(Int.(readdlm("rslts/aims_sim_n$(prefixn)_$(prefixp)_$(k)pops_path_$(s)"))[:])
    #    SnpArrays.filter("n$(prefixn)_$(prefixp)_$(k)pops_ordered", trues(n), aims; des = "rslts/sim_path_n$(prefixn)_$(prefixp)_$(k)pops_path_$(s)")
    #end
end

configs = [#(1000, 10000, 5, "1k", "10k"),
    #(1000, 100000, 5, "1k", "100k"),
    (1000, 1000000, 5, "1k", "1m"),
    #(10000, 10000, 5, "10k", "10k"),
    #(10000, 100000, 5, "10k", "100k"),
    #(1000, 10000, 10, "1k", "10k"),
    #(1000, 100000, 10, "1k", "100k"),
    #(1000, 1000000, 10, "1k", "1m"),
    #(10000, 10000, 10, "10k", "10k"),
    #(10000, 100000, 10, "10k", "100k"),
]
for config in configs
    println(config)
    run_case(config...)
    println("done $config")
end
