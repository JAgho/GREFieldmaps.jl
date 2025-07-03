using ArgParse, SpheriCart, Krylov, NIfTI, StatsBase, Statistics, LinearAlgebra, StaticArrays, SphericalHarmonicExpansions, MultivariatePolynomials


function nifti_coordinates(nifti::NIVolume)
    dims = size(nifti)
    aff = getaffine(nifti)
    carts = SVector{4, Float64}.(Tuple.(CartesianIndices((dims... ,1))))
    coords = SVector{3, Float64}.(map(x->(aff*x)[1:3] , carts))[:,:,:,1]
    fmap_mid_svec = SVector{3, Float64}((aff * SVector{4, Float64}(((dims .÷ 2) .+ 1)..., 1))[1:3])
    svec_coords = [SVector{3, Float64}(x) .- fmap_mid_svec for x in coords] ./1000
    return svec_coords
end

function map_to_tesla(hzmap, b0=63.3e-3, gamma = 42.58e6)
    f0 = b0 * gamma
    return (hzmap .+ f0) ./ gamma
end
function compute_smooth_map(fname)
    L = 5
    b0 = niread(fname)
    mask = .!(isnan.(b0))

    linds = findall(vec(mask))
        
    svec_coords = nifti_coordinates(b0)
    fmap_vals = b0[linds]

    basis = SolidHarmonics(L)
    sph, gr = compute_with_gradients(basis, vec(svec_coords))
    A = sph[linds, :]
    b = map_to_tesla(fmap_vals)
    (outcome, stats) = Krylov.lsqr(A, b)
    recon_field = reshape(sph * outcome, size(b0))
    #recon_field[linds] .= b
    return recon_field
end

function compute_smooth_coeffs(fname)
    L = 5
    b0 = niread(fname)
    mask = .!(isnan.(b0))

    linds = findall(vec(mask))
        
    svec_coords = nifti_coordinates(b0)
    fmap_vals = b0[linds]

    basis = SolidHarmonics(L)
    sph, gr = compute_with_gradients(basis, vec(svec_coords))
    A = sph[linds, :]
    b = map_to_tesla(fmap_vals)
    (outcome, stats) = Krylov.lsqr(A, b)
    #recon_field = reshape(sph * outcome, size(b0))
    #recon_field[linds] .= b
    return outcome, sph
end

function coeff_to_grad(outcome, coords)
    @polyvar x y z
    c = SphericalHarmonicCoefficients(outcome, 1.0, false)
    f = sphericalHarmonicsExpansion(c, x, y, z)
    sph = @fastfunc f
    
    fx = differentiate(f, x)
    fy = differentiate(f, y)
    fz = differentiate(f, z)
    gx = @fastfunc fx
    gy = @fastfunc fy
    gz = @fastfunc fz
    
    grad_field = map(x -> SVector{3, Float64}(
    Base.invokelatest(gx, x...), 
    Base.invokelatest(gy, x...), 
    Base.invokelatest(gz, x...)
	), coords)
    return grad_field
end

function build_A_matrix(coords, mask, L, gradscheme::AbstractMatrix{Float})

end



mask = .!(isnan.(niread("data/b0_nofac_1ms.nii")))
antimask = .!(mask)
b0 = compute_smooth_map("data/b0_nofac_1ms.nii")
xg = compute_smooth_map("data/xshim_nofac_redo_1ms.nii")
yg = compute_smooth_map("data/yshim_nofac_1ms.nii")
zg = compute_smooth_map("data/zshim_nofac_1ms.nii")
xgrad = (xg .- b0)
xgrad[.!(mask)] .= NaN
ygrad = (yg .- b0)[:,24,:]
ygrad[.!(mask)] .= NaN
zgrad = (zg .- b0)
zgrad[.!(mask)] .= NaN

contour(xgrad[:,24,:], aspect_ratio=1, levels=50, title = "x-Gradient")
contour(ygrad[24,:,:], aspect_ratio=1, levels=50, title = "y-Gradient")
contour(zgrad[:,24,:], aspect_ratio=1, levels=50, title = "z-Gradient")

plot_3plane(xgrad)
coords = nifti_coordinates(niread("data/b0_nofac_1ms.nii"))
xgrad_linear  = getindex.(coords, 3).* (0.2e-3)
contour(xgrad_linear[:,24,:], aspect_ratio=1, levels=50, title = "x-Gradient Linear")
contour((xgrad .- xgrad_linear)[:,24,:].*1e5, aspect_ratio=1, levels=50, title = "x-Gradient Nonlinear")
contour((xgrad)[:,:,24].*1e5, aspect_ratio=1, levels=50, title = "x-Gradient Nonlinear")


b0_c, sph = compute_smooth_coeffs("data/b0_nofac_1ms.nii")
x_c, sph_x = compute_smooth_coeffs("data/xshim_nofac_redo_1ms.nii")

x_c .- b0_c

recon_x = reshape(sph * (x_c .- b0_c), size(mask))

recon_x[.!(mask)] .= NaN
heatmap(recon_x[:,24,:].*1e3, aspect_ratio=1, levels=50, title = "x-Gradient Linear")

#contour(diff(recon_x, dims=1)


grad_field = coeff_to_grad(x_c.-b0_c, coords)
#test =  @SVector [NaN, NaN, NaN]
xtemp = getindex.(grad_field, 1)
xtemp[antimask] .= NaN
analytic = xtemp*1e5
finite = diff(((xgrad )),dims=1).*(1e5*(1/0.0046)) 
heatmap(analytic[:,24,:], aspect_ratio=1, levels=50, title = "x-Gradient finite")

heatmap(finite.- maximum(finite), aspect_ratio=1, levels=50, title = "x-Gradient analytic")

