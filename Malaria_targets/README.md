# Hands-on session: APO-HOLO systems for malaria targeting
## Introduction and context
This documents will provide a brief overview of the APO (without ligand) and HOLO (with ligand) protein kinase G systems that is given as exercise for the exam and how to set up the simulation box to get the simulations running.

Protein kinase G (PKG) is a key regulatory enzyme found in the malaria parasite Plasmodium. In Plasmodium falciparum, the species responsible for the most severe form of malaria in humans, PKG plays a critical role in various stages of the parasite's lifecycle, including its development within the mosquito vector and its replication within the human host's red blood cells. PKG is involved in signalling pathways that regulate essential processes such as parasite invasion of host cells, egress from infected cells, and differentiation into sexual forms for transmission to mosquitoes. Additionally, PKG has been identified as a potential target for antimalarial drug development due to its essential role in parasite survival and its divergence from human homologs, making it an attractive candidate for selective inhibition. The corresponding UNIPROT is [A5K0N4](https://www.uniprot.org/uniprotkb/A5K0N4/entry). The PDB IDs are [5DYL](https://www.rcsb.org/structure/5DYL) for the APO form, [5DZC](https://www.rcsb.org/structure/5DZC) for the HOLO form bound to a nonhydrolyzable analog of adenosine triphosphate (called ANP in the simulation box), [5F0A](https://www.rcsb.org/structure/5F0A) for the HOLO form bound to an inhibitor (called 1FB in the simulation box), and [5FET](https://www.rcsb.org/structure/5FET) for the HOLO form bound to another inhibitor (called 1TR in the simulation box).

The exact role of these molecules when targeting PKG is not strictly relevant for this tutorial. If you are interested you can read the corresponding publication[^1] from which the APO and HOLO + ANP structures were taken, while the other two come from results still under peer revision and, therefore, unpublished. Each excercise directory contains two subdirectories, one with the APO form (which is the same for all exercises) and one with one of the three ligand-bound structures linked above. These are not membrane proteins, and as such they can be simulated in a box of pure water, analogously to what done during the Lysozyme tutorial.

The individual directories of the different protein/ligand systems are organised as in the following
```
PKG_APO/HOLO_X_excercise
│   ionize.mdp
│   reference_topology_PKG_XXX.gro
│   reference_topology_PKG_XXX.top
│   sbatch_me.sh
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
- The `ionize.mdp`, `reference_topology_PKG_XXX.gro`, and `reference_topology_PKG_XXX.top` files. These are the building blocks to set up your starting configuration. These will be covered in depth later in the tutorial.
- A `sbatch_me.sh` file. This is a text file that contains the set of instructions that will run your starting configuration through energy minimisation, NVT, NPT, and eventually the production phase.
- A `forcefield` directory. Inside this, you can find several text files that define 'numerically' the key components of your system. You might recognise some of them, like `PKG.itp`, the water model `tip3.itp`, and the ligands `1FB.itp`, `1TR.itp`, and `ANP.itp`. You can peek inside them and take a look at how a molecule is described numerically in MD simulations. The format of these files - that is, how these numbers are organised and reported - depends on the software you are using. Since in this course you use GROMACS, these files are written in 'GROMACS' format. However, while these numbers and these files might be different for other softwares, the core idea of numerically representating with specific functions the intra- and inter-molecular interactions is paramount in molecular dynamics.
- A `solutions_files` directory. Here there are the `start.gro` and `topol.top` files that you should obtain at the end of this tutorial. If you want to run the simulations instantly, then move the files from this directory directly to `PKG_APO/HOLO_X_excercise` and go at the last sections of this tutorial about running the simulation.
- The `step1_em`, `step2_nvt`, and `step3_npt` directories. Inside each of these you can find one or more `mdp` (molecular dynamics parameters) files that will be used to compile and run the energy minimization, the NVT, and the NPT equilibrations. These will be run sequentially by the script `sbatch_me.sh`, taking your starting configuration through a set of equilibration phases to relax the starting configuration gradually and avoid the system exploding.
- The `step4_prod` directory. Here there is the `prod.mdp` file, which is the final mdp file for the production run of the system. This is usually the longest part of the simulation, which can take 3 to 5 days on Baobab depending on the system. The final output files of the simulation will be inside this directory, and the analysis for the exam will mainly revolve around the `prod.xtc` trajectory file.

## A look into the starting files and the general aim of this guide

The idea of this tutorial is to give you the main ingredients to build your own solvated simulation box without going through the hassle of finding a good target to simulate, fix the experimental structure files, find a consistent force field to describe the system, and embedding the protein in the membrane. Considering this, inside the main directory you can find a copy of three files which will be the basic blocks to build the system, namely
- `ionize.mdp`, the mdp file to add the ions in the system.
- `reference_topology_PKG_XXX.gro`, the starting structure that contains the structure of PKG in its APO form or in its HOLO form alongside the corresponding docked ligand.
- `reference_topology_PKG_XXX.top`, the starting topology of the system.

From these, and by using as reference the [Lysozyme in water](http://www.mdtutorials.com/gmx/lysozyme/) tutorial for solvation and ion addition, you should be able to obtain two files, the starting configuration of the solvated protein, which you will call `start.gro`, and the updated topology of the system, which you will call `topol.top`. These files will be then picked up by the `sbatch_me.sh` script and used as starting point to run the equilibration and the production runs of the system.

In the following, the tutorial shows the set-up of the simulation box by using as example the HOLO+ANP structure. Nevertheless, the logic of the procedure is identical and can be applied directly to all the other configurations as well.

## Building the starting box
First of all, send the directory with the exercises to your home in Baobab (with `scp`) and login into Baobab. Then, request a node for a couple of hours with `salloc` in the following way
```
salloc --ntasks=1 --cpus-per-task=4 --partition=private-gervasio-cpu --time=120:00
```
Finally, source the GROMACS installation
```
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load GROMACS/2023.1-CUDA-11.7.0
```
and verify that the sourcing was okay by typing
```
gmx --version
```
Notice how here you are requesting 4 CPUs but no GPU, differently from one of the first lectures (in fact, you are on `--partition=private-gervasio-cpu`). In this case, you are going to use the allocation on Baobab just to set up the box and not to run the simulation, so the GPU is not needed. Once the allocation has been granted to you, change directory and go into one of the two exam systems by using `cd`. Here, if you list what's inside (type `ls`), you should see the various `stepX_XXX` directories and the starting files.

You are now ready to assemble the starting configuration of your system. One last point before moving on. GROMACS has a lot (> 100) tools that are accessible by typing `gmx` followed by a keyword. In this tutorial you will use `solvate`, `grompp`, `genion`, and `make_ndx`. If you have any doubts remember that you can look online for the explanation of the tool and which flags are needed (`-f`, `-o`, `-s` etc.). For example, [this](https://manual.gromacs.org/current/onlinehelp/gmx-solvate.html) is the manual page of `gmx solvate`. All this information is also awailable on the spot if you type `-h` or `--h` (for *help*) after the tool's name, e.g., `gmx solvate --h`.

### Solvating and generating the box
At the beginning, you can take a look at the starting topology (`reference_topology_PKG_HOLO_ANP.top`) and the starting configuration to solvate (`reference_topology_PKG_HOLO_ANP.gro`). The topology reads like this
```
#include "./forcefield/forcefield.itp"
#include "./forcefield/PKG.itp"
#include "./forcefield/ANP.itp"
#include "./forcefield/tip3.itp"
#include "./forcefield/ions.itp"

[ system ]
PKG HOLO with ANP and MG

[ molecules ]
; Compound         #mols
PKG                    1
MG                     1
ANP                    1
SOL                    2
```
First of all, notice that sometimes there are lines starting with a semi-column `;`. These are **comments**, that is, lines that are not read by GROMACS. Comments are used to annotate the files and write down details that are helping you - the users - remember what you are doing, what is the meaning of some variables, etc. To all effects, these are equivalent to your notes on the border of a book or on a slide to write something that is worth remembering, e.g. some explanation of the professor. They are ignored by the software and should be informative for the person writing or for the person supposed to read the code. Feel free to write your own comments, if it helps you. Just remember to put a `;` at the beginning of each line that you do not want to be read by GROMACS (you can't forget, as if you do the GROMACS command will simply fail and complain about non-sensical commands).

Following up, the `#include` statements tell GROMACS where to find the elements of the force field. As can be seen, they are collected inside the `forcefield` directory and must appear with a specific order. First, the set of parameters defining the force field (`forcefield.itp`). Then, the definition of the individual molecules (`PKG.itp`, `ANP.itp`, etc.) that will populate you system. These do not have to appear in a specific order, however <ins>all</ins> the molecules that you intend to use <ins>must</ins> be defined here. Building the topology, that is, filling in this file, is roughly the equivalent of running `gmx pdb2gmx` on the pdb file of the protein, as you did in the Lysozyme tutorial. In this case the system is much more complex and the pdb needs further handling. Thus, it wouldn't be possible to do in one simple line within GROMACS, so the topology is already provided. After the `#include` statements, there is the `[ system ]` section, which is simply the name of the system. It is worth giving the system a meaningful name to help in recognising the systems in the future.

Finally, there is the `[ molecules ]` section. This is a very important section which must contain <ins>all</ins> the molecules of the system <ins>in the order in which they appear</ins>. Notice that, for the time being, the topology contains one protein (`PKG`), one magnesium ion (`MG`), one ligand (`ANP`), and two water molecules (`SOL`). This is unique for the HOLO system with ANP, because the ANP ligand is an adenosine triphosphate analogue which is coordinated to a magnesium ion together with two water molecules. This is important because the whole complex ANP + MG + 2 H<sub>2</sub>O is the correct physiological state of the HOLO state of this protein when bound to ANP. The ligand doesn't bind without the coordinating water molecules and magnesium ion. GROMACS doesn't know if a ligand necessitas a magnesium or other ions and water molecules in some given positions. Thus, they are part of the initial box **before** the solvation as they have been modelled already in the correct state. This is where the biological and chemical knowledge of the user is critical for building a physiologically correct system. Among the systems presented in this tutorial, ANP is the only ligand that necessitates other molecules to be modelled properly in its binding site. The topologies and starting structures of the other systems contain only the host protein PKG and a ligand (for the HOLO states) or the protein alone (for the APO case). The APO structure and the three ligands are shown in Figure 1.

| ![Figure 1](../images/apoandligands.png) |
|:--:|
| Figure 1 *From left to right: PKG in its APO form; binding site with ANP, the magnesium ion, and the two coordinating water molecules; binding site with 1FB; binding site with 1TR. The residues nearest to the ligands are represented as lines. Color legends: red (oxygen), blue (nitrogen), yellow (phosphorus), pink (magnesium). Carbons are colored in cyan for the protein and iceblue for the ligands. Hydrogen atoms are hidden exception made for the water molecules, where they are white.* |

Now, take a look at the starting configuration `reference_topology_PKG_HOLO_ANP.gro` by opening this file with a text reader. The `.gro` file has a fixed format, and it is better to not modify it by hand if you are not completely sure about what you are doing. The first lines look like this
```
PKG HOLO with ANP ligand
12908
    1ACE     H1    1  12.255   9.683   1.920
    1ACE    CH3    2  12.282   9.770   1.861
    1ACE     H2    3  12.301   9.855   1.927
    1ACE     H3    4  12.202   9.795   1.792
    1ACE      C    5  12.405   9.741   1.784
    1ACE      O    6  12.461   9.632   1.791
    2GLY      N    7  12.452   9.840   1.705
    2GLY      H    8  12.403   9.929   1.702
[...]
```
The file has a first line which contains the title of the box (`PKG HOLO with ANP ligand`), a second line which contains the number of the atoms in the box (`12908`), and then it contains in order all the atoms of the system. These are organised usually as nine columns. The first (here `1ACE`) is the specific number and name of the reidue - which in this case is the capped N-terminal of PKG. The second column contains the specific name of the atom, and usually the first letter indicates the element (here you have a hydrogen followed by a carbon, two hydrogen atoms, another carbon, and so on). The third is simply the number of the entry. It always starts with `1` and goes up to the number of elements in the box. Then, columns four to six contain the x, y, and z coordinates of that atom, while the columns seven to nine contain its velocity, reported by axial component. Here you can see how the velocities are not reported, and so you have only six columns. All atoms must always have a position, but might have zero or undefined velocity. This is the case now, as you are building the box from a static experimental image. One of the main roles of the equilibrations phase is this - to relax the starting positions and assign reasonable starting velocities to all the atoms.

At the other end of the file, the last lines look like this
```
[...]
  800ANP    HC812901   9.954  11.756   3.126
  800ANP   HN3B12902  10.448  11.664   2.794
  801SOL     OW12903  10.681  11.334   3.167
  801SOL    HW112904  10.690  11.235   3.182
  801SOL    HW212905  10.712  11.356   3.073
  802SOL     OW12906  10.542  11.544   3.341
  802SOL    HW112907  10.543  11.515   3.439
  802SOL    HW212908  10.540  11.644   3.339
   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
```
The meaning of the columns is the same as before. The last molecules to appear are the two water molecules (and in fact, in the topology, they are reported last). The last line of a `.gro` file has, like the first two, a special meaning. It contains the coordinates of the box, that is, the length of the box along x, y, and z. Sometimes it can have more than three numbers for particular box shapes. You can see how this box is actually undefined - GROMACS keeps this line used, as it has a special role, but it is filled with zeros. You will need to define it before solvating the system.

As in the Lysozyme tutorial, you can set up the box by using `gmx editconf`
```
gmx editconf -f reference_topology_PKG_HOLO_ANP.gro -o reference_topology_PKG_HOLO_ANP_boxed.gro -bt dodecahedron -c -d 1
```
Here you are sending it the reference structure and asking GROMACS to position your protein in the centre (`-c`), built a dodecahedron box (`-bt dodecahedron`), and make sure that the protein is at least 1nm distant (`-d 1`) from the box sides to avoid periodic boundary conditions (PBCs) artifacts. Please not that, exactly because of PBCs, there is not a real centre of the box, and positioning the protein in it with `-c` is mostly for representation purposes rather than real physics. 

Differently from the Lysozyme tutorial, here you are also specifying the shape type of the box. In terms of simulations results, the outcome should not depend on the dimension or shape of the simulation box, as long as the box is large enough to accomodate the contents (the protein in this case) and avoid self-interaction between PBC images. So, why bother changing the shape or dimension? Why don't just make a huge cube? A hint should come from one of the last lines of the ouput of `gmx editconf`, which will look something like this
```
new box volume  :1496.67               (nm^3)
```
Try now to build the same system by specifying a cubic (`-bt cubic`) box, you will end up with something like this
```
new box volume  :2116.61               (nm^3)
```
This means that you will have a roughly ~40% larger box, which means more water to solvate, a higher number of molecules to simulate, and therefore slower simulations. An example of how different box types look like is reported in Figure 2.

| ![Figure 2](../images/boxtypes.png) |
|:--:|
| Figure 2 *Box types* |

Before running `gmx solvate`, you have to know which water model you want to use. For this force field, as also reported in the `reference_topology_PKG_HOLO_ANP.top` file, the model is [TIP3P](https://en.wikipedia.org/wiki/Water_model) (you are importing the parameters with this line `#include "./forcefield/tip3.itp"`), a three-point water model. This means that each water molecule in your simulation will simply be represented with three sites: one oxygen atom and two hydrogen atoms. To access a three-points water model, the flag name for `gmx solvate` is `-cs spc216.gro`, the same as for the Lysozyme tutorial.

Thus, you can solvate the system with the following
```
gmx solvate -cp reference_topology_PKG_HOLO_ANP_boxed.gro -cs spc216.gro -o reference_topology_PKG_HOLO_ANP_solvated.gro
```
where you are asking to add water to the structure with `-cp reference_topology_PKG_HOLO_ANP_boxed.gro`, use as reference a three-points water model with `-cs spc216.gro`, add call the resulting output structure `-o reference_topology_PKG_HOLO_ANP_solvated.gro`. Some of you may notice that in the Lysozyme tutorial you also have to pass the topology of the system with the `-p` flag. This is not mandatory, and I prefer not to do it because the updated topology will have the same name as the input one, which is prone to errors. However, you will have to update by hand the topology to include the presence of the water molecules.

The (last lines of the) output of this command will look something like this
```
[...]
Volume                 :     1496.67 (nm^3)
Density                :     999.599 (g/l)
Number of solvent molecules:  44578 
```
GROMACS tries to fill the box with water to reach the density of ca. 1g/L. The most important part here is the number of water molecules inserted, in this case `44578` (this number can oscillate slightly). You now has just to add this number at the end of the topopolgy to tell GROMACS that the content of the box has changed.

First, copy the reference topology and call it `reference_topology_PKG_HOLO_ANP_solvated.top` with the `cp` (copy) command
```
cp reference_topology_PKG_HOLO_ANP.top reference_topology_PKG_HOLO_ANP_solvated.top
```
Then, add the amount of water molecules to `reference_topology_PKG_HOLO_ANP_solvated.top` by correcting the `[ molecules ]` section in the following way
```
[...]
[ molecules ]
; Compound         #mols
PKG                    1
MG                     1
ANP                    1
SOL                44580
```
You have to add the name of the water in the box (`SOL`) and the number reported after cleaning up the system with the python script. For the HOLO + ANP structure, `SOL` was already reported as there are already two water molecules coordinating the MG ions, and as such you just have to add two to the total solvation number. For the other structure, just report the plain number obtained from the output of `gmx solvate`.

At this point, the structure of the solvated system is contained in `reference_topology_PKG_HOLO_ANP_solvated.gro`, while its topology in `reference_topology_PKG_HOLO_ANP_solvated.top`.

## References
[^1]: El Bakkouri, M., et al. "Structures of the cGMP-dependent protein kinase in malaria parasites reveal a unique structural relay mechanism for activation." Proceedings of the National Academy of Sciences 116.28 (2019): 14164-14173. [DOI:10.1073/pnas.190555811](https://doi.org/10.1073/pnas.190555811)

| ![Alt Text](https://media.giphy.com/media/vFKqnCdLPNOKc/giphy.gif) |
|:--:|
| Figure 2 *gif example* |
