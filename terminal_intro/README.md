# **A brief introduction to terminals**

## General considerations

Let's start with a couple of definitions.

- **Terminal** A terminal is a program that provides a text-based interface to your computer. Trivially, it is the window where you type commands and see the output, acting as the communication interface between you and the system. The standard in macOS is [Terminal.app](https://en.wikipedia.org/wiki/MacOS)), in Linux Ubuntu it is the [GNOME Terminal](https://en.wikipedia.org/wiki/GNOME_Terminal), while in Windows it is generally the [PowerShell](https://en.wikipedia.org/wiki/PowerShell).
- **Shell** A shell is the command interpreter that runs inside the terminal. It reads the commands you type, executes them, and displays the results. Common examples include [bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) and [zsh](https://en.wikipedia.org/wiki/Z_shell).

You can see the relationship between a human user, a terminal, a shell, and a computer as two people talking in a room. One person is you, the other is the computer, and the room where you are talking is the terminal—a common space where you can directly interact. The shell is the common set of rules you both understand: your shared language. Depending on the operating system and your preferences, you can “speak” with a computer in different rooms and in different languages. Generally, *terminal* and *shell* are used interchangeably, although they are not the same.

Before the advent of the most common graphical user interfaces, most of the interactions with computers happened thorough a terminal with a keyboard. Nowadays, instead, most people interact with computers using a combination of keyboard and mouse/touch-pad by clicking on windows, files, etc. This is equivalent to using a terminal, simply the commands are interpreted as clicks or specific keystrokes, e.g., double-clicking on a directory to open it or pressing `Del` to delete a file.

In scientific sectors, the usage of the terminal is still very common. There are several reasons for this. One is that many scientific tools have been developed for being used within terminals. Another is that basically all super computers (like Baobab from the University of Geneva) are usually accessed online and have no graphical interface, only a terminal one. Most importantly, using a terminal (especially in a [Unix](https://en.wikipedia.org/wiki/Unix)/[Unix-like](https://en.wikipedia.org/wiki/Unix-like) environment like macOS or Linux Ubuntu) provides a universal common ground that works the same on laptops, servers, and high‑performance computing systems.

This tutorial is a brief introduction to common terminal vocabulary, that is, we will learn **commands** for a given **shell**. Here, we focus on commands available in both `bash` and `zsh`, two extremely widespread and popular shells.

---

# Part 1 — Getting Comfortable

Start by opening your terminal apps.  
Since this guide focuses on `bash` and `zsh`:

- macOS users → open **Terminal.app**  
- Linux users → open the terminal distributed with your OS  
- Windows users → **do not** use PowerShell (different language). Instead open the Windows Subsystem for Linux ([WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux)).

If we set-up together the WSL app, you should see that you are using some version of Linux Ubuntu as the first message when you launch it. Formally, for you this is equivalent as being on another operative systems (Linux Ubuntu instead of Windows) and you are already inside the terminal, as no graphical interface is shipped with WSL.

As few general notes before starting with the tutorial.

 - When you read that you have to `run` a command, it means that you have to type it within the terminal and to press `Enter`.
 - Remember that you can always check online how to use a command and, by running `command --help` you can access a brief explanation of the tool (e.g. `cp --help`).
 - You can use the `Tab` key on your keyboard for fast autocompletion of paths and commands.

---

## **Exercise 1: Where Am I?**  
**Goal:** Learn to check the current directory and it's content.

**Commands:** `pwd`, `ls`

Without a graphical interface, you cannot see where you are located (Desktop? Downloads? Home? etc.) and what is inside the directory. You can check these by using the `pwd` (*print working directory*) and `ls` (*list*) commands, respectively.

a) Run `pwd`  

The command will answer with something of the form like the following
```bash
/directory1/directory2/final_directory
```
This means that, in this moment, you are inside `final_directory` (the last entry of the line), which itself is inside `directory2` that is inside `directory1`. Notice how the names of the different directories are separated by the symbol `/`. The core idea is that there is a starting directory, called **root** (indicated by the first symbol `/` on the left) that containts all directories. Then, like a tree from its root, all the following content of the computer is inside directories of directories. This d is the standard organization of all computers (altough in Windows Terminals like Powershell you might see `\` rather than `/` to separate directories).

Very importantly, the output of `pwd` is the **path** to your directory, that is, a unique combination of the position and name that gives the location of the directory in the computer.

b) Run `ls`

This command will show the content of the directory where you are (that is, the one given by `pwd`) as a list of files. If it shows nothing, then the directory is empty.

---

## **Exercise 2: Navigation**  
**Goal:** Learn how to change directories, and absolute & relative paths.

**Commands:** `cd`

To change a directory, you run `cd` (*change directory*) and insert the path of the directory you want to move to. Generally, when you open a terminal, you are in your `home` directory. From wherever you are, you can always run

a) `cd ~`, or

b) `cd`

to get back to your home. Try it, you should not be moving from where you are.
Now, let's suppose you have a directory here called `Documents` (it should be visible with `ls`). To access it you run

c) `cd Documents`

or, more generally, you can access a directory that is in your current directory with the command

d) `cd name_directory`

To move one directory up, i.e., go into the directory that contains the one where we are now, you can use

e)  `cd ..`  

where the `..` is a general way to say *one up*. You can concatenate the points, i.e., `cd ../..` means go **two** directories up, `cd ../../..` means **three** and so on. The single dot `.`, insted, indicates *this* directory where you are now. For example, running

f) `cd .`

moves you...nowhere, since you are telling the computer to move to where you are to where you are. Similarly,

g) `ls .`

will list the files and directories present in the current working directory, equivalently to `ls`. While the single dot might seem (and it is) useless in the context of these commands, it might be useful for others that you will encounter down the road.

This is a good moment to grasp the concept of paths. As said before, a **path** is a text string that tells the computer where a file or directory is located in the filesystem. Example of paths are

`/home/user/Documents`
`Documents`
`./notes/todo.txt`
`../data/results.csv`

Path can be  **absolute** or **relative**. An **absolute path** *always* starts from the root of the filesystem, represented by `/`. It is the full address of a file or directory and is valid no matter where you currently are. A **relative path**, instead, is defined in relation to your current directory, and it *never* starts with `/`.

In the four examples above, the first is an absolute path to the `Documents` directory, while all the others are relative paths. This means that, wherever you are, if you type `cd /home/user/Documents` you will *always* be taken to that specific `Documents` directory. Typing `cd Documents` instead will take you inside the directory `Documents` contained in your current working directory, or will fail if the directory is not present. Similarly, `./notes/todo.txt` and `../data/results.csv` are relative paths to files, that is, the first means *I am referring to the `todo.txt` file inside the directory `notes` inside `.` (this) directory* and the second *I am referring to the `results.csv` file inside the `data` directory contained in `..` (one directory above me)*.

To make the difference more clear and state better the importance of this difference, let's suppose you have two different `Documents` directories with absolute paths `/home/bob/Documents` and `/home/alice/Documents`. The command `cd Documents` will work both in `bob` and `alice` directories, but will take you to two directories that have the same name but *different* absolute paths, that is, they are *different directories* that will have different files, sub-directories etc. You can easily imagine that absolute paths are **unique**, as you can't have files or directories with same name within the same folders (as in all computers).

Finally,

h) Try navigating using both absolute and relative path. Remember that you can always go back to your home with `cd` or `cd ~`.

---

# Part 2 — Creating Things

## **Exercise 3: Create directories**  
**Goal:** Learn how to create a directory

**Commands:** `mkdir`

To create a new directory, you can run the command `mkdir` (*make directory*). Specifically, running `mkdir name_directory` will generate a directory with name `name_directory` in the directory where you are, equivalently to `mkdir ./name_directory`. You can use asbsolute and relative paths as well, e.g., `mkdir ../name_directory` to create it inside the directory containing the current one or `mkdir /home/bob/Documents/name_directory` to create it specifically at that path.

a) Create a directory called `playground`

b) Enter this directory

c) Inside it, create two directories: `animals` and `plants`. You can create multiple directories at once with `mkdir name_dir1 name_dir2 name_dir3 etc...`.

---

## **Exercise 4: Make Some Files**  
**Goal:** Learn how to generate files 
**Commands:** `touch`

You can create files with the `touch` command. Enter the `animals` directory:

a) Create three files: `dog.txt`, `cat.txt`, `elephant.txt`

b) Check that they exist using `ls`

c) What is their relative path? And their absolute one?

d) Go back two directories to the one containing `playground` with `cd ../..` and try to create `bamboo.txt` and `baobab.txt` inside `plants` **without** entering in `plants` but by using relative or absolute paths, something like `touch [path_to_plants]/bamboo.txt`

---

# Part 3 — Viewing & Editing

## **Exercise 5: Read Your Files**
**Goal:** Learn how to visualize a file

**Commands:** `cat`, `echo`

You can access the content of a file by using the command `cat`. The simple version `cat file` will print the content of `file` on the screen. `echo`, instead, is used for printing something specific. `echo "Hello World!"` will print `Hello World!` on the screen.

a) Navigate to the directory containing the file `dog.txt`

b) Display the contents of `dog.txt` using `cat` run `cat dog.txt`. It should be empty.

c) Try the command `echo "Hello World!"`

d) Add a line using to the file using `echo "Dogs are friendly animals." > dog.txt`

You will notice that the `echo` command is unchanged, but we added a `> dog.txt` at the end of it. This means *whatever is the output of the command before `>`, put it in the file `dog.txt`.

---

## **Exercise 6: Append Content**  
**Goal:** Learn output redirection (>>)

a) Append another sentence to `dog.txt` using `>>`, e.g. `echo "your cool sentence" >> dog.txt`

b) Display the updated file with cat.

Notice that `>` **rewrites** the content of the file, that is, whatever is contained in the files gets overwritten with the new input, while `>>` **appends** to the content of the file, i.e., it adds the input at the end of it.

---

# Part 4 — Moving and Renaming Things

## **Exercise 7: Rename a File**
**Goal:** Learn how to rename and move files and directories

**Commands:** `mv`

You can rename or move a file with the command `mv` (*move*). Example usages are `mv a.txt b.txt`, which renames file `a.txt` in `b.txt`, and `mv a.txt ./Documents/`, which moves file `a.txt` from the current directory inside the `Documents` directory. You can also combine them like `mv a.txt ./Documents/b.txt` which both moves and renames the file. The same can be applied to directories.

a) Rename the `elephant.txt` into `toad.txt` in the animal folder

b) Create `otter.dat` in the directory `plants` and move it in the directory `animals`

c) Rename the `plants` directory in `animal_garden` and move it inside the `animals` directory. Try to do it in one command.

---

# Part 5 — Copying and Removing Things

## **Exercise 8: Make a Copy**
**Goal:** Learn how to make a copy of a file or a directory

**Commands:** `cp`

You can copy files with the command `cp` (*copy*), e.g. `cp a b` will copy fle `a` in `b`. Remember that absolute and relative paths can be used everywhere we specify a file, for example `cp ../a.txt /home/Documents/a_backup.txt` will copy the `a.txt` file from the directory one above us to a file named `a_backup.txt` inside the directory `/home/Documents/`

a) Copy the `dog.txt` file both into `dog.txt` and `dog_backup.txt`. Which one works? Why?

b) Copy the full directory `animals` with `cp animals animals_2`. Does it work?

As you will notice, directories can't be simply copied as they might contain other files, and you have to specify what to do with them. To copy a full directory you have to specify the flag `-r` (for *recursive*) so that the copy command knows that the order has to be repeated for everything contained inside the directory as well.

c) Copy the full directory `animals` with `cp -r animals animals_2`. Does it work? What's inside the directory `animals_2`?

## **Exercise 9: Delete Things**
**Goal:** Learn how to delete a file or a directory

**Commands:** `rm`, `rmdir`

You can remove files and directories with `rm` (*remove*) and `rmdir` (*remove directory*). **IMPORTANT: There is no Recycling Bin. Everything that is deleted is lost, permanently. Double think before deleting stuff, and double check your commands. Removing files and directories from certain System directories can permanently damage the operative system of your computer, or you can easily delete one month worth of data analysis. As a safety rule, do not touch things that you are not sure what they are.**

a) Create a directory called `poaching` inside `playground`

b) Remove it with `rmdir`, e.g. `rmdir poaching`

c) Create it again, and inside it create two files, `rifle.txt` and `poacher.txt`

d) Remove the directory with `rmdir`. Does it work? Why?

e) Now remove it with `rm`, e.g. `rm poaching`. Does it work? Why?

f) Remove the `rifle.txt` file with the `rm` command

g) And now remove the full `poaching` directory with the `rm -r` command. Here the `-r` flag option of `rm` has the same value as in `cp`, i.e., *recursive*, so apply this command to this directory and recursively to everything that is inside it, that is, remove this directory and everything inside.

# Part 6 — Text editors
**Goal:** A first approach to text editors in a terminal

**Commands:** `vim`, `nano`

Until now, you saw how to add data to a file by using the `echo` command. This is a neat trick to store on file output from commands, but in general it is not very usefull to open and modify files. There are several text editors available in terminals, and which one to use largely depends on the personal preferences. Two of the most common are [`vim`](https://www.vim.org/) and [`nano`](https://en.wikipedia.org/wiki/GNU_nano). These are complex software packages with a billion options, and sadly the learning curve is a bit steep, but the basics for opening, modifying and saving (a basic usage of a text editor) are very simple.

a) Open the `toad.txt` file in the `animal` folder with `vim` by running `vim toad.txt`

b) Go into insert mode by pressing `a` (you will see the last line of the page change into **--INSERT--**)

c) Type whatever you want. Remember, you have to move with arrows and keyoard as clicking on the terminal is not available

d) Press `Esc`. The **--INSERT--** mode is not active anymore. You can now access the general options by typing `:` (you will see the symbol appear at the bottom of the page, outside the file text)

e) To save the file, write `:w` (*write*) and press enter. Now you can quit with `:q` (*quit*).

f) Re-open the file with `vim`. Go again in `INSERT` mode and type something. Now try to close the file *without* saving with `:q`. Does it work? Why?

g) You can force the closure with the `!`symbol, e.g., `:q!` will close without saving any changes.
