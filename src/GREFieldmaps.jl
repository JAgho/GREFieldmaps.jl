module GREFieldmaps
    println("loading environment...")
using ArgParse, SpheriCart, Krylov, NIfTI, StatsBase, Statistics, LinearAlgebra, StaticArrays, SphericalHarmonicExpansions, MultivariatePolynomials

export process_nifti, main

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

function process_nifti(fmap_path::String, target_path::String, basename::String, L::Int)
    fmap = niread(fmap_path)
    target = niread(target_path)
    println("reading input NIfTIs...")
    mask = fmap .> 0
    linds = findall(vec(mask))
    
    svec_coords = nifti_coordinates(fmap)
    target_coords = nifti_coordinates(target)
    mask_coords = svec_coords[linds]
    fmap_vals = fmap[linds]
    
    basis = SolidHarmonics(L)
    sph, gr = compute_with_gradients(basis, vec(svec_coords))
    A = sph[linds, :]
    b = map_to_tesla(fmap_vals)
    println("solving inverse problem...")
    (outcome, stats) = Krylov.lsqr(A, b)
    recon_field = reshape(sph * outcome, size(fmap))
    recon_field[linds] .= b
    println("solving forwards problem...")
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
	), target_coords) .* 1000

    b0_field = map(x -> Base.invokelatest(sph, x...), target_coords) .* 1000
    println("saving NIfTIs...")
    b0_nii = NIVolume(target.header, target.extensions, b0_field)
    grad_x = NIVolume(target.header, target.extensions, getindex.(grad_field,1))
    grad_y = NIVolume(target.header, target.extensions, getindex.(grad_field,2))
    grad_z = NIVolume(target.header, target.extensions, getindex.(grad_field,3))
    
    niwrite("$(basename)_b0.nii", b0_nii)
    niwrite("$(basename)_grad_x.nii", grad_x)
    niwrite("$(basename)_grad_y.nii", grad_y)
    niwrite("$(basename)_grad_z.nii", grad_z)
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--fmap"
            help = "Path to the fieldmap NIfTI file"
            arg_type = String
            required = true
        "--target"
            help = "Path to the target NIfTI file"
            arg_type = String
            required = true
        "--basename"
            help = "Base name for output NIfTI volumes"
            arg_type = String
            required = true
	"--L"
            help = "Spherical harmonic order L to model up to"
            arg_type = Int
            default = 7
    end
    
    args = parse_args(ARGS, s)
    process_nifti(args["fmap"], args["target"], args["basename"], args["L"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


end # module GREFieldmaps
