# papago
Patch based gaussian mixture modelling of medical imaging  


## Pre-processing
- First, we need data to be aligned, and we need the masks of those alignments (perhaps even the interp-matrices?)

## Execution

#### Data Preparation
Preprocess subvolumes, e.g. `adniprep.m`, get `md` struct. This involves building the `medialDataset` structure via `restorationmd` and processing the images via `processmd`.

#### Training On Cluster
1. Break up dataset into "subvolume columns" via grid, e.g. via `md2subvols.m`, and save subvolume columns. This can take a long amount of time and/or memory.
1. Run `wgmm` on cluster distributed on each subvolumes, e.g. via `sgeTrain.sh`.
This can be done via `model0` (isotropic data) of `model3` (weighted data).

#### Testing On Cluster
1. run `mccRecon.m` (`papago.recon`) on all patches on each subvolume, and store the reconstructions!
1. (unfinished) re-compose volume.

#### Evaluation for optimal K at each location (On Cluster)
Loop steps for Training and Testing for various K. Choose K based on best patch reconstruction at each location.

#### Training and Testing On a Single Machine
This is usually done on a subset of the image grid.  
Loop over (sub)grid:
1. create/load subvolume column.
1. run `wgmm` via `papago.train`
1. run `papago.recon` to reconstruct all the patches in this subvolume (perhaps for just a subset of subjects)  

Quilt patches.
