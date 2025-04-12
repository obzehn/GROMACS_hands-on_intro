# Hands-on session: G protein-coupled receptors (GPCRs)
## Introduction and context
This documents will provide a brief overview of the [GPCRs](https://en.wikipedia.org/wiki/G_protein-coupled_receptor) systems given as exercise and how to set up the simulation box to get the simulations starting.

The excercise directory contains two subdirectories, each one with a slightly different structure of the same transmembrane domain of an adhesion GPCR[^1]. These are membrane proteins, and as such your starting configuration consists in the GPCR embedded in a mixed lipid bilayer with cholesterol and phospholipids (POPC).

The individual directories of the different structures are organised as in the following

```
GPCR_structure_X_excercise
│   ionize.mdp
│   reference_topology_GPCR_structure_1.gro
│   reference_topology_GPCR_structure_1.top
│   sbatch_me.sh
│   water_removal.py
│
└─── forcefield
│   │   [various forcefields files]
│
└─── solution_files
│   │   [the solutions of the excercise]
│
└─── step1_em
│   │   em.mdp
│
└─── step2_nvt
│   │   [various nvt.mdp files]
│
└─── step3_npt
│   │   [various npt.mdp files]
│
└─── step4_prod
    │   prod.mdp
```
Going by order, you can see
- The `ionize.mdp`, `reference_topology_GPCR_structure_1.gro`, `reference_topology_GPCR_structure_1.top`, and `water_removal.py` files. These are the building blocks to set up your starting configuration. We will cover these in depth later in this tutorial.
- A `sbatch_me.sh` file. This is a text file that contains the set of instructions that will run your starting configuration through energy minimisation, NVT, NPT, and eventually the production phase.
- A `forcefield` directory. Inside this, you can find several text files that define 'numerically' the key components of your system. You might recognise some of them, like `GPCR.itp`, `cholesterol.itp`, and `POPC.itp`. You can peek inside them and take a look at how a molecule is described numerically in MD simulations. The format of these files - that is, how these numbers are organised and reported - depends on the software you are using. Since in this course we use GROMACS, these files are written in 'GROMACS' format. However, while these numbers and these files might be different for other softwares, the core idea of numerically representating with specific functions the intra- and inter-molecular interactions is paramount in molecular dynamics.
- A `solutions_files` directory. Here there are the `index.ndx`, `start.gro`, and `topol.top` files that you should obtain at the end of the tutorial. If you want to run the simulations instantly, then move the files from this directory directly to `GPCR_structure_X_excercise` and go at the last section of this tutorial about running the simulation.
- The `step1_em`, `step2_nvt`, and `step3_npt` directories. Inside each of these you can find one or more `mdp` (molecular dynamics parameters) files that will be used to compile and run the energy minimization, the NVT, and the NPT equilibrations. These will be run sequentially by the script `sbatch_me.sh`, taking your starting configuration through a set of equilibration phases to relax the starting configuration gradually and avoid the system exploding.
- The `step4_prod` directory. Here there is the `prod.mdp` file, which is the final mdp file for the production run of the system. This is usually the longest part of the simulation, which can take 3 to 5 days on Baobab depending on the system. The final output files of the simulation will be inside this directory, and the analysis for the exam will mainly revolve around the `prod.xtc` trajectory file.
