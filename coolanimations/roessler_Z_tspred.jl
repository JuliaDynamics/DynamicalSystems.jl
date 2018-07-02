using Plots
using TimeseriesPrediction
using DynamicalSystemsBase
#pyplot()
#gr()

cd(@__DIR__)

ds = Systems.roessler()
tr = trajectory(ds, 2000)

start_pred = 50000
strain = tr[1:start_pred,3]
spred = localmodel_tsp(strain, 50, 15, 7000)



@time anim = @animate for i = start_pred-500:25:start_pred+6000
    println(i, " / ", (start_pred-500:25:start_pred+6000)[end])
    l = @layout [a b]

    # plot gray background Attractor
    p1 = path3d(columns(tr)... , xlim=(-10,12), ylim=(-11,8), zlim=(0,23),
        #xlab = "x", ylab = "y", zlab = "z",
        title = "Roessler Attractor",
        seriesalpha = 0.1,
        seriescolor = :black,
        label = "")

    # draw red dot with tail on top
    path3d!(p1, columns(Dataset(tr[i-800:i]))...,
            seriescolor= :black, seriesalpha =0.5, lw=1, label="")
    path3d!(p1, [tr[i,1]], [tr[i,2]], [tr[i,3]], seriescolor= :red, label="",
    marker=1, markerstrokecolor= :red)

    # right plot
    p2 = plot(tr[i-1600:i+400, 3],
        xlim=(0,2000), ylim=(0,23), label="Z Variable",
        seriescolor=:green,
        xticks=nothing,
        xlabel="Time",
        ylabel="Z Coordinate")
    scatter!(p2, [1600], [tr[i, 3]], label="", lw = 3,seriescolor=:green, markerstrokecolor= :green)
    #vline!(p2, [800], line = (:black, 1), label="")

    j = i - start_pred
    vline!(p2, [1600 - j], line = (:black, 1), label="")

    #add predicted stuff
    if j > 0
        plot!(p2, (1601-j):2000,spred[1:j+400], label="Prediction",   seriescolor=:red)
        scatter!(p2, [1600],   [spred[j]], label="", lw = 3,seriescolor=:red, markerstrokecolor= :red)
        scatter!(p2, [1601-j], [spred[1]], label="", lw = 3,seriescolor=:red, markerstrokecolor= :red)
    end
    p = plot(p1, p2, size = (600,300), layout= l)

end
#gif(anim,randstring() * ".gif")

mp4(anim, "roessler_Z_tspred.mp4")
