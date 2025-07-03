midslice(x::AbstractVector{<:Number}, matrix_size::Tuple{Int, Int, Int}; dims=2) = 
    midslice(reshape(x, matrix_size...), dims=dims)

function plot_3plane(a, matrix_size::Tuple{Int, Int, Int})
    x = reverse(midslice(a, matrix_size, dims=1), dims=1)
    y = midslice(a, matrix_size, dims=2)
    z= reverse(midslice(a, matrix_size, dims=3)')
    plt1 = heatmap(x, aspect_ratio=1, color=:grays, cbar=false, title="Coronal",axis=([], false))
    plt2 = heatmap(y, aspect_ratio=1, color=:grays, cbar=false, title="Axial",axis=([], false),)
    plt3 = heatmap(z, aspect_ratio=1, color=:grays, cbar=false, title="Sagittal",axis=([], false),)
    plt = plot(plt1, plt2, plt3, layout=(1, 3), size=(1200, 400))
    display(plt)
end

function plot_3plane(a)
    x = reverse(midslice(a, dims=1), dims=1)
    y = midslice(a, dims=2)
    z= reverse(midslice(a, dims=3)')
    plt1 = heatmap(x, aspect_ratio=1, color=:grays, cbar=false, title="Coronal",axis=([], false))
    plt2 = heatmap(y, aspect_ratio=1, color=:grays, cbar=false, title="Axial",axis=([], false),)
    plt3 = heatmap(z, aspect_ratio=1, color=:grays, cbar=false, title="Sagittal",axis=([], false),)
    plt = plot(plt1, plt2, plt3, layout=(1, 3), size=(1200, 400))
    display(plt)
end