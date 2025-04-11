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
- The `step1_em`, `step2_nvt`, and `step3_npt` directories. Inside each of these you can find one or more `mdp` (molecular dynamics parameters) files that will be used to compile and run the energy minimization, the NVT, and the NPT equilibrations. These will be run sequentially by the script `sbatch_me.sh`, taking your starting configuration through a set of equilibration phases to relax the starting configuration gradually and avoid the system exploding.
- The `step4_prod` directory. Here there is the `prod.mdp` file, which is the final mdp file for the production run of the system. This is usually the longest part of the simulation, which can take 3 to 5 days on Baobab depending on the system. The final output files of the simulation will be inside this directory, and the analysis for the exam will mainly revolve around the `prod.xtc` trajectory file.

## A look into the starting files and the general aim of this guide

The idea of this tutorial is to give you the main ingredients to build your own solvated simulation box without going through the hassle of finding a good target to simulate, fix the experimental structure files, find a consistent force field to describe the system, and embedding the protein in the membrane. Considering this, inside the main directory you can find a copy of four files which will be the basic blocks to build the system, namely
- `ionize.mdp`, the mdp file to add the ions in the system.
- `reference_topology_GPCR_structure_1.gro`, the starting structure that contains the trans-membrane domain of the GPCR already embedded in the membrane.
- `reference_topology_GPCR_structure_1.top`, the starting topology of the system.
- `water_removal.py`, a python script that removes the water molecules that GROMACS wrongly positions within the lipid bilayer during the solvation phase.

From these, and by using as reference the [Lysozyme in water](http://www.mdtutorials.com/gmx/lysozyme/) tutorial for solvation and ion addition, you should be able to obtain three files, the starting configuration of the solvated membrane, which we will call `start.gro`, the updated topology of the system, which we will call `topol.top`, and the index file of the system, which we will call `index.ndx` and that is useful for controlling complex multi-phase systems like a protein embedded in a membrane bilayer. These files will be then picked up by the `sbatch_me.sh` script and used as starting point to run the equilibration and the production runs of the system. This tutorial uses the first GPCR structure as an example, but the process is identical for all the four configurations.

## Building the starting box
First of all, send the directory with the exercises to your home in Baobab (with `scp`) and login into Baobab. Then, request a node for a couple of hours with `salloc` in the following way
```
salloc --ntasks=1 --cpus-per-task=4 --partition=private-gervasio-gpu --time=120:00
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
Notice how here we are requesting 4 CPUs but no GPU, differently from one of the first lectures (in fact, we are on `--partition=private-gervasio-cpu`). In this case, we are going to use the allocation on Baobab just to set up the box and not to run the simulation, so the GPU is not needed. Once the allocation has been granted to you, change directory and go into one of the two exam systems by using `cd`. Here, if you list what's inside (type `ls`), you should see the various `stepX_XXX` directories and the starting files.

You are now ready to assemble the starting configuration of your system. One last point before moving on. GROMACS has a lot (> 100) tools that are accessible by typing `gmx` followed by a keyword. In this tutorial you will use `solvate`, `grompp`, `genion`, and `make_ndx`. If you have any doubts remember that you can look only for the explanation of the tool and which flags are needed (`-f`, `-o`, `-s` etc.). For example, [this](https://manual.gromacs.org/current/onlinehelp/gmx-solvate.html) is the manual page of `gmx solvate`. All this information is also awailable on the spot if you type `-h` or `--h` (for *help*) after the tool's name, e.g., `gmx solvate --h`.

### Solvating and cleaning up residual water
At the beginning, we can take a look at the starting topology (`reference_topology_GPCR_structure_1.top`) and the starting configuration to solvate (`reference_topology_GPCR_structure_1.gro`). The topology reads like this
```
#include "./forcefield/forcefield.itp"
#include "./forcefield/cholesterol.itp"
#include "./forcefield/GPCR.itp"
#include "./forcefield/POPC.itp"
#include "./forcefield/tip4pd.itp"
#include "./forcefield/ions.itp" 

[ system ]
; Name
GPCR structure 1

[ molecules ]
; Compound        #mols
GPCR                1
CHL    	           70
POPC  	          280
```
First of all, notice that sometimes there are lines starting with a semi-column `;`. These are **comments**, that is, lines that are not read by GROMACS. We use these to annotate the files and write down comments that are helping us - the users - remember what we are doing, what is the meaning of some variables, etc. To all effects, these are equivalent to your notes on the border of a book or on a slide to write something that is worth remembering, e.g. some explanation of the professor. They are ignored by the software and should be informative for the person writing or for the person supposed to read the code. Feel free to write your own comments, if it helps you. Just remember to put a `;` at the beginning of each line that you do not want to be read by GROMACS (you can't forget, as if you do the GROMACS command will simply fail and complain about non-sensical commands).

Following up, the `#include` statements tell GROMACS where to find the elements of the force field. As can be seen, they are collected inside the `forcefield` directory and must appear with a specific order. First, the set of parameters defining the force field (`forcefield.itp`). Then, the definition of the individual molecules (`cholesterol.itp`, `GPCR.itp`, etc.) that will populate our system. These do not have to appear in a specific order, however <ins>all</ins> the molecules that we intend to use <ins>must</ins> be defined here. Building the topology, that is, filling in this file, is roughly the equivalent of running `gmx pdb2gmx` on the pdb file of the protein, sa you did in the Lysozyme tutorial. In this case the system is much more complex and it wouldn't be possible to do in one simple line within GROMACS, so the topology is already provided. After the the `#include` statements, there is the `[ system ]` section, which is simply the name of the system. It is worth giving the system a meaningful name to help in recognising the systems in the future.

Finally, there is the `[ molecules ]` section. This is a very important section which must contain <ins>all</ins> the molecules of the system <ins>in the order in which they appear</ins>. Notice that, for the time being, the topology contains one GPCR, 70 cholesterol molecules (`CHL`), and 280 phospholipids (`POPC`).

Now, take a look at the starting configuration `reference_topology_GPCR_structure_1.gro` by opening this file with a text reader. The `.gro` file has a fixed format, and it is better to not modify it by hand if you are not completely sure about what you are doing. The first lines look like this
```
GPCR structure 1
47140
    1NTHR     N    1   5.174   5.111   3.908  0.0000  0.0000  0.0000
    1NTHR    H1    2   5.156   5.041   3.838  0.0000  0.0000  0.0000
    1NTHR    H2    3   5.086   5.138   3.950  0.0000  0.0000  0.0000
    1NTHR    H3    4   5.223   5.065   3.984  0.0000  0.0000  0.0000
    1NTHR    CA    5   5.256   5.221   3.857  0.0000  0.0000  0.0000
[...]
```
The file has a first line which contains the title of the box (`GPCR structure 1`), a second line which contains the number of the atoms in the box (`47140`), and then it contains in order the all the atoms of the system. These are organised usually as nine columns. The first (here `1NTHR`) is the specific number and name of the reidue - which in this case is the N-terminal of the GPCR. The second column contains the specific name of the atom, and usually the first letter indicates the element (here we have a nitrogen followed by three hydrogen and a carbon and so on). The third is simply the number of the entry. It always starts with `1` and goes up to the number of elements in the box. Then, columns four to six contain the x, y, and z coordinates of that atom, while the columns seven to nine contain its velocity, reported by axial component. Notice how all atoms always have a position, but might have zero velocity. This is out case now, as we are building the box from a static experimental image. One of the main roles of the equilibrations phase is this - to relax the starting positions and assign reasonable starting velocities to all the atoms.

At the other end of the file, the last lines look like this
```
[...]
 1193POPC  C11847137   9.898   1.890   4.009  0.0000  0.0000  0.0000
 1193POPC  H18R47138   9.885   1.793   3.961  0.0000  0.0000  0.0000
 1193POPC  H18S47139   9.807   1.907   4.069  0.0000  0.0000  0.0000
 1193POPC  H18T47140   9.985   1.888   4.075  0.0000  0.0000  0.0000
  10.06535  10.06535  10.56089
```
The meaning of the columns is the same as before. The last molecules to appear are the phopsholipids (and in fact, in the topology, they are reported last). The last line of a `.gro` file has, like the first two, a special meaning. It contains the coordinates of the box, that is, the length of the box along x, y, and z. For example, this box is roughly a cube, with x and y lengths of 10.06535 nm and z length of 10.56089 nm. Sometimes, for special types of boxes (like those used for the protein-ligand simulations), there are more that three values reported. The starting configuration, which I am visualising with [VMD](https://www.ks.uiuc.edu/Research/vmd/), is reported in Figure 1.

| ![Figure 1](../images/startgpcr.png) |
|:--:|
| Figure 1 *Side (left) and top (right) view of the system. The GPCR is reported in red as cartoons, while the lipids are shown in sticks colored in red (oxygen), blue (nitrogen), cyan (carbon), white (hydrogen), and yellow (phosphorus). The box is shown as blue lines.* |

The first step is then to solvate the system. In the Lysozyme tutorial, before solvation, the `gmx editconf` command is used to prepare the box and to make it large enough to fit the protein. Here, we do not use it as the box is already prepared and it depends on the width of the lipid bilayer.

Before running `gmx solvate`, we have to know which water model we want to use. For this force field, as also reported in the `reference_topology_GPCR_structure_1.top` file, the model is [TIP4PD](https://en.wikipedia.org/wiki/Water_model) (we import the parameters with this line `#include "./forcefield/tip4pd.itp"`), a four-point water model. This means that each water molecule in your simulation will be actually represented with four sites: one oxygen atom, two hydrogen atoms, and a dummy site near the oxygen atom that has negative charge. This dummy atom is a numerical trick to better represent the distribution of charge of the lone pair in water's oxygen. To access four-points water models, the flag name for `gmx solvate` is `-cs tip4p.gro`.

Thus, you can solvate the system with the following
```
gmx solvate -cp reference_topology_GPCR_structure_1.gro -cs tip4p.gro -o GPCR_structure_1_solvated.gro
```
where you are asking to add water to the structure with `-cp reference_topology_GPCR_structure_1.gro`, use as reference a four-points water model with `-cs tip4p.gro`, add call the resulting output structure `-o GPCR_structure_1_solvated.gro`. Some of yoy may notice that in the Lysozyme tutorial you also have to pass the topology of the system with the `-p` flag. This is not mandatory, and I prefer not to do it because the updated topology will have the same name as the input one, which is prone to errors. However, you will have to update by hand the topology to include the presence of the water molecules.

The (last lines of the) output of this command will look something like this
```
Volume                 :     1069.94 (nm^3)
Density                :        1013 (g/l)
Number of solvent molecules:  21053
```
GROMACS tries to fill the box with water to reach the density of ca. 1g/L. The most important part here is the number of water molecules inserted, in my case `21053` (this number can oscillate slightly). In principle, you could add this number at the end of the topopolgy and you would be done with the solvation, as in the Lysozyme tutorial. However, you might think that the system has actually more than one phase (lipids and protein) and that GROMACS tried to insert water molecules in eaxh and every empty volume it could find. How did it behave near the protein? And how near the lipids? You can check with VMD.

The results for my box are shown in Figure 2. Water has been in all places where it was possible, also between the lipids, which we know is not realistic as this is a strong lipophilic (non-hydrophobic) region. Generally, this is not a major problem, and within the first nanoseconds of run the water molecules get expelled naturally because of the chemical nature of the lipid tail region. Nevertheless, for delicate systems like this that involve proteins highly susceptible to not physiological starting configurations, it is better to get rid of them and clean the system before adding the ions.

You can remove these water molecules with the `water_removal.py` Python script, which you can run with the following command
```
python3 water_removal.py GPCR_structure_1_solvated.gro GPCR_structure_1_solvated_clean.gro
```
The script takes as input `GPCR_structure_1_solvated.gro`, check where there are overlapping water molecules with the lipids and the protein, removes them, and returns the cleaned system under the name `GPCR_structure_1_solvated_clean.gro`. For me, the output looks like this
```
Initial number of water molecules: 21053
Number of water molecules to be deleted: 2117
Final number of water molecules: 18936
```
The script removed ca. 2000 water molecules. You can see now, in Figure 2, that the system has been cleaned and looks much more biologically sound.
| ![Figure 2](../images/solvatedgpcr.png) |
|:--:|
| Figure 2 *On the left, the solvated system after `gmx solvate`. On the right, the same system after removing ca. 2000 water molecules with the `water_removal.py`. Protein and lipids are reported in grey, while water molecules in red (oxygen) and white (hydrogen). The fourth site of this water model is in a lighter red tone, but it nearly overlaps with the water oxygen atoms, and is thus nearly invisible.* |


topology update

image before, after solvation, after cleaning

### Adding ions
gmx genion --maxwarn 1

topology update

### Generate the index file

[^1]: I. Liebscher, T. Schöneberg, and D. Thor. "Stachel-mediated activation of adhesion G protein-coupled receptors: insights from cryo-EM studies." Signal transduction and targeted therapy 7.1 (2022): 227. [DOI:10.1038/s41392-022-01083-y](https://doi.org/10.1038/s41392-022-01083-y)
