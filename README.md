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

Papers
--------
If you find this library useful, please cite ([download bib](citations.bib)):  

- **Medical Image Imputation from Image Collections**  
[A.V. Dalca](http://adalca.mit.edu), [K.L. Bouman](https://people.csail.mit.edu/klbouman/), [W.T. Freeman](https://billf.mit.edu/), [M.R. Sabuncu](http://sabuncu.engineering.cornell.edu/), [N.S. Rost](https://www.massgeneral.org/doctors/doctor.aspx?id=17477), [P. Golland](https://people.csail.mit.edu/polina/)  
IEEE TMI: Transactions on Medical Imaging 38.2 (2019): 504-514. eprint [arXiv:1808.05732](https://arxiv.org/abs/1808.05732)  

- **Population Based Image Imputation**  
[A.V. Dalca](http://adalca.mit.edu), [K.L. Bouman](https://people.csail.mit.edu/klbouman/), [W.T. Freeman](https://billf.mit.edu/), [M.R. Sabuncu](http://sabuncu.engineering.cornell.edu/), [N.S. Rost](https://www.massgeneral.org/doctors/doctor.aspx?id=17477), [P. Golland](https://people.csail.mit.edu/polina/)  
In Proc. IPMI: International Conference on Information Processing and Medical Imaging. LNCS 10265, pp 1-13. 2017. 
