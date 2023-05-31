Dataset = DynamicalSystems.Dataset

plot_dataset(args...; kwargs...) = plot_dataset!(Scene(), args...; kwargs...)
function plot_dataset!(scene, data::Dataset{2}, color = :black; kwargs...)
    makiedata = [Point2f(d) for d in data]
    scatter!(scene, makiedata; color = color, markersize = 0.01, kwargs...)
    return scene
end
function plot_dataset!(scene, data::Matrix, color = :black; kwargs...)
    size(data, 2) != 2 && error("Matrix dimension 2 must be 2")
    makiedata = [Point2f(data[i, 1], data[i, 2]) for i in 1:size(data, 1)]
    scatter!(scene, makiedata; color = color, markersize = 0.01, kwargs...)
    return scene
end
function plot_dataset!(scene, data::Dataset{3}, color = :black; kwargs...)
    makiedata = [Point3f(d) for d in data]
    lines!(scene, makiedata; color = color, transparency = true, linewidth = 2.0, kwargs...)
    return scene
end

plot_datasets(args...; kwargs...) = plot_datasets!(Scene(), args...; kwargs...)
function plot_datasets!(scene, datasets, colors; kwargs...)
    length(datasets) == length(colors) || error("Datasets and colors must "*
    "have equal length")

    for i in eachindex(datasets)
        plot_dataset!(scene, datasets[i], colors[i]; kwargs...)
    end
    return scene
end
