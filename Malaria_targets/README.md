# Hands-on session: APO-HOLO systems for malaria targeting
## Introduction and context
This documents will provide a brief overview of the APO (without ligand) and HOLO (with ligand) protein kinase G systems that is given as exercise for the exam and how to set up the simulation box to get the simulations running.

Protein kinase G (PKG) is a key regulatory enzyme found in the malaria parasite Plasmodium. In Plasmodium falciparum, the species responsible for the most severe form of malaria in humans, PKG plays a critical role in various stages of the parasite's lifecycle, including its development within the mosquito vector and its replication within the human host's red blood cells. PKG is involved in signalling pathways that regulate essential processes such as parasite invasion of host cells, egress from infected cells, and differentiation into sexual forms for transmission to mosquitoes. Additionally, PKG has been identified as a potential target for antimalarial drug development due to its essential role in parasite survival and its divergence from human homologs, making it an attractive candidate for selective inhibition. The corresponding UNIPROT is [A5K0N4](https://www.uniprot.org/uniprotkb/A5K0N4/entry). The PDB IDs are [5DYL](https://www.rcsb.org/structure/5DYL) for the APO form, [5DZC](https://www.rcsb.org/structure/5DZC) for the HOLO form bound to a nonhydrolyzable analog of adenosine triphosphate (called ANP in the simulation box), [5F0A](https://www.rcsb.org/structure/5F0A) for the HOLO form bound to an inhibitor (called 1FB in the simulation box), and [5FET](https://www.rcsb.org/structure/5FET) for the HOLO form bound to another inhibitor (called 1TR in the simulation box).

The exact role of these molecules when targeting PKG is not strictly relevant for this tutorial. If you are interested you can read the corresponding publication[^1] from which the APO and HOLO + ANP structures were taken, while the other two come from results still under peer revision and, therefore, unpublished. Each excercise directory contains two subdirectories, one with the APO form (which is the same for all exercises) and one with one the three ligand-bound structures linked above. These are not membrane proteins, and as such they can be simulated in a box of pure water, analogously to what done during the Lysozyme tutorial.

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

## References
[^1]: El Bakkouri, M., et al. "Structures of the cGMP-dependent protein kinase in malaria parasites reveal a unique structural relay mechanism for activation." Proceedings of the National Academy of Sciences 116.28 (2019): 14164-14173. [DOI:10.1073/pnas.190555811](https://doi.org/10.1073/pnas.190555811)
