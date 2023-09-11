using Distributions
using OpenADMIXTURE, HaploADMIXTURE
using DataFrames, CSV
using StatsBase, Random
using SnpArrays
using StableRNGs
using Mmap
using DelimitedFiles

using LinearAlgebra

LinearAlgebra.BLAS.set_num_threads(1)

function draw_allele_freqs(i, K, freqdat, fst; rng=Random.GLOBAL_RNG)
    freqs = freqdat[:,i]
    fst_val1, fst_val2 = fst[2(i-1)+1], fst[2i]
    
    fst = 2 * fst_val1 * fst_val2 / (fst_val1 + fst_val2)
    r = rand(rng, Dirichlet((1 - fst) / fst .* [freqs[1], freqs[2], freqs[3], freqs[4]]), K)
    HaploADMIXTURE.project_p!(r', [1, 2, 3, 4], K; pseudocount = 0.005)
    return r'
end

# n_region: number of regions
# K: number of populations
# n_ind: number of individuals to sample
# n_samples_per_region: number of samples per region
# alpha: dispersion of cluster (the lower, the further)
# gamma: dispersion within cluster (the higher, the closer)

function main()
    n_ind = parse(Int, ARGS[1])
    n_snps = parse(Int, ARGS[2])
    n_region = 50
    n_samples_per_region = n_ind ÷ n_region
    K = parse(Int, ARGS[3])
    gamma = 50.0
    chunk_size = 1000

    seed = parse(Int, ARGS[4])
    alpha = parse(Float64, ARGS[5])
    rng = StableRNG(seed)

    # sample q for each region
    region_q = rand(rng, Dirichlet(alpha .* ones(K)), n_region)
    # generate full q matrix
    q = mapreduce(i -> rand(rng, Dirichlet(gamma * region_q[:, i] .+ 0.00001), n_samples_per_region), hcat, 1:n_region)
    OpenADMIXTURE.project_q!(q, collect(1:K)) # projection to min 1e-5 simplex. 
    dat = CSV.read("/home/kose/HaploADMIXTURE/tgp_superpop_param.txt", DataFrame, delim=' ', ignorerepeated=true)

    freqdat = readdlm("/home/kose/HaploADMIXTURE/tgp_chr1_freqs.txt")
    inds = sample(rng, 1:(size(dat, 1) ÷ 2), n_snps, replace=true)
    for (i, ind) in enumerate(inds)
        while true
            fst1 = dat.FST[2(ind-1)+1]
            fst2 = dat.FST[2ind]
            if fst1 <= 0.0001 || fst1 > 0.9999 || fst2 <= 0.0001 || fst2 > 0.9999
                inds[i] = sample(rng, 1:size(dat, 1) ÷ 2)
                ind = inds[i]
            else
                break
            end
        end
    end

    p = mapreduce(j -> draw_allele_freqs(inds[j], K, freqdat, dat.FST; rng=rng), hcat, 1:(n_snps ÷ 2))


    writedlm("q_$(K)pops_$(n_ind)indiv_$(n_snps)snps_$(seed)_$(alpha).txt", q)
    writedlm("p_$(K)pops_$(n_ind)indiv_$(n_snps)snps_$(seed)_$(alpha).txt", p)


    s = SnpArray("g_$(K)pops_$(n_ind)indiv_$(n_snps)snps_$(seed)_$(alpha).bed", n_ind, n_snps)
    nchunks = n_snps ÷ 2chunk_size
    for c in 1:nchunks
        local_freq = q' * p[:, (c - 1) * chunk_size * 4 + 1: c * chunk_size * 4]
        for i in 1:n_ind
            for j in 1:chunk_size
                mnd = Multinomial(2, (local_freq[i, 4(j-1) + 1: 4(j-1) + 4]))
                genotype = rand(rng, mnd)
                genotype1 = genotype[1] + genotype[2]
                genotype2 = genotype[1] + genotype[3]
                s[i, (c - 1) * 2chunk_size + 2(j-1)+1] = (genotype1 == 0 ? 0x00 : genotype1 == 1 ? 0x02 : 0x03)
                s[i, (c - 1) * 2chunk_size + 2(j-1)+2] = (genotype2 == 0 ? 0x00 : genotype2 == 1 ? 0x02 : 0x03)
            end
        end
    end

    Mmap.sync!(s.data)

    open("g_$(K)pops_$(n_ind)indiv_$(n_snps)snps_$(seed)_$(alpha).fam", "w") do io
        for i in 1:n_ind
            write(io, join([i, i, 0, 0, 0, 0], "\t") * "\n")
        end
    end

    open("g_$(K)pops_$(n_ind)indiv_$(n_snps)snps_$(seed)_$(alpha).bim", "w") do io
        for j in 1:n_snps
            write(io, join([1, j, 0, j, 1, 2], "\t") * "\n")
        end
    end
end

main()
