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

As a general note, when you read that you have to `run` a command, it means that you have to type it within the terminal and to press `Enter`.

---

## **Exercise 1: Where Am I?**  
**Goal:** Learn to check the current directory and it's content.
Commands: `pwd`, `ls`

Without a graphical interface, you cannot see where you are located (Desktop? Downloads? Home? etc.) and what is inside the directory. You can check these by using the `pwd` (*print working directory*) and `ls` (*list*) commands, respectively.

a) Run `pwd`  

The command will answer with something of the form like the following
```
/directory1/directory2/final_directory
```
This means that, in this moment, you are inside `final_directory` (the last entry of the line), which itself is inside `directory2` that is inside `directory1`. This is the standard organization of all computers. The core idea is that there is a starting directory, called **root** (indicated by the first symbol `/`) that containts all directories. Then, like a tree from its root, all the following content of the computer is inside directories of directories.

b) Run `ls`

This command will show the content of the directory where you are (that is, the one given by `pwd`) as a list of files. If it shows nothing, then the directory is empty.

---

## **Exercise 2: Navigating Like a Pro**  
**Goal:** Learn `cd`, absolute & relative paths.

a) Move to your home directory: `cd ~`  
b) Navigate to the `Documents` folder  
c) Move one directory up using `cd ..`  
d) Try navigating using both an absolute and a relative path  

---

# Part 2 — Creating Things

## **Exercise 3: Create a Playground**  
**Goal:** Learn `mkdir`, `cd`.

a) Create a directory called `playground`  
b) Enter this directory  
c) Inside it, create two directories: `animals` and `plants`

---

## **Exercise 4: Make Some Files**  
**Goal:** Learn `touch`.

Inside the `animals` directory:

a) Create three files: `dog.txt`, `cat.txt`, `elephant.txt`  
b) Check that they exist using `ls`  

---

# Part 3 — Viewing & Editing

## **Exercise 5: Read Your Files**  
**Goal:** Learn `cat`.

a) Display the contents of `dog.txt` using `cat`  
b) Add a line using:  
   ```bash
   echo "Dogs are friendly animals." > dog.txt
