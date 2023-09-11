using OpenADMIXTURE, HaploADMIXTURE
using Random, StableRNGs, DelimitedFiles, SnpArrays
using CUDA
using JuMP
using HiGHS
function jac!(m, q1, q2; cutoff = 1e-5)
    I, K = size(q1)
    @assert size(m) == (K, K)
    @assert size(q2) == (I, K)
    @inbounds for k2 in 1:K
        for k1 in 1:K
            cnt = 0
            s = zero(eltype(m))
            for i in 1:I
                if q1[i, k1] > cutoff || q2[i, k2] > cutoff
                    cnt += 1
                    s += (q1[i, k1] - q2[i, k2]) ^ 2
                end
            end
            m[k1, k2] = 1 - sqrt(s / 2cnt)
        end
    end
    m
end
function jac(q1, q2; cutoff=1e-5)
    I, K = size(q1)
    m = Matrix{Float64}(undef, K, K)
    jac!(m, q1, q2; cutoff=cutoff)
end
function permute_q2(q1, q2; cutoff=1e-5)
    I, K = size(q1)    
    m = jac(q1, q2; cutoff=cutoff)
    model = Model(HiGHS.Optimizer)
    @variable(model, 0 <= x[1:K, 1:K] <= 1)
    @objective(
        model,
        Max,
        sum(m[i, j] * x[i, j] for i in 1:K, j in 1:K)
    )
    @constraint(model, sum(x[1:K, j] for j in 1:K) .== 1)
    @constraint(model, sum(x[i, 1:K] for i in 1:K) .== 1)
    optimize!(model)
    perm = value.(x)
    collect(transpose(perm * transpose(q2)))
end
rmse_(q1, q2) = sqrt(sum((q1 .- q2).^2)/prod(size(q1)))

function main()
    prefix = ARGS[1]
    qfile = ARGS[2]

    I = parse(Int, ARGS[3])
    J = parse(Int, ARGS[4]) ÷ 2
    K = parse(Int, ARGS[5])
    rng = StableRNG(3517)
    d_oa = OpenADMIXTURE.run_admixture(prefix * "_oa.bed", K; use_gpu=true, rng=rng)

    T = Float64
    T2 = Float64
    g = SnpLinAlg{T2}(SnpArray(prefix * "_ha.bed"))
    Q = 3
    d_ha = HaploADMIXTURE.AdmixData2{T, T2}(I, J, K, Q, g; rng=rng)
    d_cu, g_cu = HaploADMIXTURE._cu_admixture_base(d_ha, g, d_ha.I, d_ha.J)
    HaploADMIXTURE.init_em!(d_ha, g, 5; d_cu=d_cu, g_cu=g_cu)
    HaploADMIXTURE.admixture_qn!(d_ha, g, 300, 1e-7; d_cu=d_cu, g_cu=g_cu, mode=:ZAL)

    q = readdlm(qfile)
    q_correct = q'
    function rmse(q)
        sqrt(sum(x -> x^2, q_correct .- permute_q2(q_correct, q))/prod(size(q_correct)))
    end

    rmse_oa = rmse(d_oa[1].q')
    rmse_ha = rmse(d_ha.q')

    writedlm(prefix * "_oa_p.txt", d_oa[1].p)
    writedlm(prefix * "_oa_q.txt", d_oa[1].q)
    writedlm(prefix * "_ha_p.txt", d_ha.p)
    writedlm(prefix * "_ha_q.txt", d_ha.q)
    println("$prefix: OpenADMIXTURE RMSE $(rmse_oa) / HaploADMIXTURE RMSE $(rmse_ha)")
end 

main()

