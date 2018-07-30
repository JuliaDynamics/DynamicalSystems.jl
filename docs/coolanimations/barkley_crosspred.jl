using TimeseriesPrediction
using FileIO
function barkley_periodic_boundary(T, Nx, Ny)
    a = 0.75
    b = 0.06
    ε = 1/12
    D = 1/1.28
    u = zeros(Nx, Ny)
    v = zeros(Nx, Ny)
    U = Vector{Array{Float64,2}}()
    V = Vector{Array{Float64,2}}()

    #Initial state that creates spirals
    u[35:end,34] = 1
    u[35:end,35] = 1
    u[35:end,36] = 1
    v[35:end,37] = 1
    v[35:end,38] = 1
    v[35:end,39] = 1


    u[1:20,14] = 1
    u[1:20,15] = 1
    u[1:20,16] = 1
    v[1:20,17] = 1
    v[1:20,18] = 1
    v[1:20,19] = 1
    v[1:20,20] = 1


    u[27:36,20] = 1
    u[27:36,19] = 1
    u[27:36,18] = 1
    v[27:36,17] = 1
    v[27:36,16] = 1
    v[27:36,15] = 1

    h = 0.75 #/ sqrt(2)
    Δt = 0.1 #/ 2
    δ = 0.01
    Σ = zeros(Nx, Ny, 2)
    r = 1
    s = 2
    function F(u, uth)
        if u < uth
            u/(1-(Δt/ε)*(1-u)*(u-uth))
        else
            (u + (Δt/ε)*u*(u-uth))/(1+(Δt/ε)*u*(u-uth))
        end
    end

    for m=1:T
        for j=1:Ny, i=1:Nx
            if u[i,j] < δ
                u[i,j] = Δt/h^2 * Σ[i,j,r]
                v[i,j] = (1 - Δt)* v[i,j]
            else
                uth = (v[i,j] + b)/a
                v[i,j] = v[i,j] + Δt*(u[i,j]^3 - v[i,j])
                u[i,j] = F(u[i,j], uth) + D * Δt/h^2 *Σ[i,j,r]
                Σ[i,j,s] -= 4u[i,j]
                Σ[  mod(i-1-1,Nx)+1,j,s] += u[i,j]
                Σ[  mod(i+1-1,Nx)+1,j,s] += u[i,j]
                Σ[i,mod(j-1-1,Ny)+1,  s] += u[i,j]
                Σ[i,mod(j+1-1,Ny)+1,  s] += u[i,j]
            end
            Σ[i,j,r] = 0
        end
        r,s = s,r
        push!(U,copy(u))
        push!(V,copy(v))
    end
    return U,V
end


D      = try Meta.parse(ARGS[1]) catch 2       end
τ      = try Meta.parse(ARGS[2]) catch 1       end
B      = try Meta.parse(ARGS[3]) catch 1       end
k      = try Meta.parse(ARGS[4]) catch 1       end
Ttrain = try Meta.parse(ARGS[5]) catch 10      end
noise  = try Meta.parse(ARGS[6]) catch 0.0      end
inp    = try ARGS[7]  catch "u"     end
outp   = try ARGS[8] catch "v"     end

println("D",D)
println("tau",τ)
println("B",B)
println("k",k)
println("Ttrain",Ttrain)
println("noise", noise)
println("inp", inp)
println("outp",outp)
c = false
w = (0,0)

Nx = 50
Ny = 50
Tskip = 300
Ttest = 50
T = Tskip + Ttest + Ttrain


U,V = barkley_periodic_boundary(T, Nx, Ny)
if inp == "v" && outp == "u"
    U, V = V, U
end

Utrain = U[1 + T - Ttest - Ttrain : T - Ttest]
noise == 0 || ( Utrain += noise * ΔU² * randn.(size.(Utrain)) )
Vtrain = V[1 + T - Ttest - Ttrain : T - Ttest]
noise == 0 || ( Vtrain += noise * ΔV² * randn.(size.(Vtrain)) )
Utest  = U[1 + T - Ttest - (D-1)τ :  T]
Vtest  = V[1 + T - Ttest - (D-1)τ :  T]

Vpred = crosspred_stts(Utrain,Vtrain,Utest, D, τ, B, k;
weighting=w, boundary=c)

fname = "Bark_$(inp)->$(outp)_D$(D)_τ$(τ)_B$(B)_k$(k)_w$(w)_tr$(Ttrain)"
noise == 0 || (fname *= "_noise$(noise)")
cd(@__DIR__); mkpath(fname); cd(fname)
save( fname * ".jld", "Vpred", Vpred, "Utest", Utest, "Vtest", Vtest)


function nrMSE(upred, utest) # 1D
    m = length(upred[1])
    p = length(upred)
    mse = 0
    for n ∈ 1:p
        mse += sum(abs2.(upred[n] - utest[n])) / p / m
    end

    umean = mean([mean(utest[i]) for i∈ 1:p ])
    nrmse = 0
    for n ∈ 1:p
        nrmse += sum(abs2.(utest[n] - umean)) / p / m
    end
    nrmse = sqrt(mse / nrmse)

    return mse, nrmse
end

inp    =  "v"
outp   =  "u"
D      =  20
τ      =  2
B      =  2
k      =  2
c      = false
w = (0,0)
Ttrain = 2000
noise = 0.0

cd(@__DIR__)

using Plots

fname = "Bark_$(inp)->$(outp)_D$(D)_τ$(τ)_B$(B)_k$(k)_w$(w)_tr$(Ttrain)"
#fname *= "_pca"
noise == 0 || (fname *= "_noise$(noise)")
cd(@__DIR__); mkpath(fname); cd(fname)


Upred, Utest, Vtest = load( fname * ".jld","Vpred", "Vtest", "Utest")
Ttest = length(Upred)
err = [abs.(Utest[i+(D-1)τ]-Upred[i]) for i=1:Ttest]

println(nrMSE(Upred,Utest[1+(D-1)τ:end]))


for i=1:10:length(Upred)
    l = @layout([a b; c d])
    p1 = plot(Vtest[i+(D-1)τ],aspect_ratio=1,st=:heatmap,seriescolor=:viridis)
    #, clims=(0,0.75)
    plot!(title = "Barkley Model")
    p2 = plot(Utest[i+(D-1)τ],aspect_ratio=1,st=:heatmap,seriescolor=:viridis)
    #, clims=(0,0.75)
    title!("original $(outp)")
    p3 = plot(Upred[i],aspect_ratio=1,st=:heatmap,seriescolor=:viridis)
    #, clims=(0,0.75),
    title!("Cross-Pred $(outp)")
    p4 = plot(err[i],aspect_ratio=1, st=:heatmap,seriescolor=:viridis)
    #clims=(0,0.1),
    title!("Absolute Error")

    p = plot(p1,p2,p3,p4, layout=l, size=(600,600))
    savefig(fname *"_$i.png")
end
