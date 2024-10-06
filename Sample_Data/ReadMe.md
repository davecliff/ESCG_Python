This folder contains some sample data to illustrate the output files written by ESCG_RPSLS.py

Be aware: the ESCG code can create a large number of output files, many of which can grow quite large.

In my implementation, the top-level data folder for an experiment is named to indicate whether it is OES or RES, how many species (NSnn), how many ablatiomns (NAn), and what the square 2D lattice side-length L is. 

So, for example, if the top-level data-folder is OES_NS05_NA0_L200 that tells me it's data for an experiment using OES with five species, no ablations, and a lattice of sidelength 200. 

Open that folder, and it has one sub-folder per computer (real or virtual) used to generate the results. If I'm using real computers, they're multi-core Apple Macs -- typically with 8, 10, or 12 CPU cores. If I'm using virtual computers, they're some kind of multi-core AWS instance -- most often I use an 128-core AWS instance, to get >100 IID results generated in parallel. 

Open the folder for any one computer, and you'll see a set of sub-folders, one per core. On a 12-core machine they might be named Core01 to Core12

Open the folder for any one core, and you'll see a set of sub-folders, one per run: a run is a specific launch of the ESCG/RPSLS simulator. The folder name for a specific launch is built from the various hyperparameter values used for that run. 

Open the folder for any one run, and you'll see the set of output files from that run: the run is a sweep of M-values so the filenames include M values in their header. The filename header can be quite long because there are several paramaters to record. For any one M-value, three CSV-format output files are written: [file_header]_domnets.csv; [file_header]_grids.csv; and [file_header]_densities.csv...

++ [file_header]_domnets.csv is the pre-ablation and post-ablation dominance networks, as binary adjacency matrices. This is a few hundred bytes.

++ [file_header]_grids.csv is a sequence of timestamped snapshots of the lattice, showing the contents of every cell in the lattice at the time of the snapshot. The frequency of the snapshots is set in the Python code, so if the *grids.csv files are growing too large for your liking, reduce the frequency of printing the snapshots. 

++ [file_header]_densities.csv has one row per Monte Carlo Step (MCS) of the simulation. On each row, the first column is the MCS number, then the number of surviving (non-extinct) species, then the number of empty cells, and then the sequence of densities for each species in the simulation.
