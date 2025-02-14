# GREFieldmaps.jl

## Overview
`GREFieldmaps` is a lightweight fitting utility for computing $B_0$ fieldmaps based on GRE-derived phase difference frequency maps. It fits frequency observations to a solid harmonic representation up to a specified basis and outputs the fitted field to conform to a specified NIfTI used as a header template. The middle of the NIfTI image is currently treated as the image isocentre.
## Installation
First install Julia, then install this package using the package mode (press `]` after starting Julia to enter package mode), and enter:

```julia
add https://github.com/JAgho/GREFieldmaps.jl
```

## Usage
Run the script using the command line:

```sh
julia --project=path/to/GREFieldmaps.jl path/to/GREFieldmaps.jl/src/GREFieldmaps.jl \
    --fmap path/to/fieldmap.nii \
    --target path/to/target.nii \
    --basename path/to/output/basename
    --L 7
```

### Arguments:
- `--fmap` : Path to the fieldmap NIfTI file.
- `--target` : Path to the target NIfTI file.
- `--basename` : Base name for output NIfTI volumes.
- `--L` : Spherical harmonic order to fit to.

This generates four output files:
- `basename_b0.nii` : Reconstructed B0 field.
- `basename_grad_x.nii` : Gradient field in the x-direction.
- `basename_grad_y.nii` : Gradient field in the y-direction.
- `basename_grad_z.nii` : Gradient field in the z-direction.

## Module Structure
- **`nifti_coordinates(nifti::NIVolume)`**: Computes real-world coordinates of a nifti file's voxels.
- **`map_to_tesla(hzmap)`**: Converts frequency map to Tesla units.
- **`process_nifti(fmap_path, target_path, basename)`**: Main processing function.
- **`main()`**: CLI interface using ArgParse.
