# **A brief introduction to terminals**

## General considerations

Let's start with a couple of definitions.

- **Terminal**  
  A terminal is a program that provides a text-based interface to your computer. Trivially, it is the window where you type commands and see the output, acting as the communication interface between you and the system.  
  The standard in macOS is [Terminal.app](https://en.wikipedia.org/wiki/MacOS)),  
  in Linux Ubuntu it is the [GNOME Terminal](https://en.wikipedia.org/wiki/MacOS),  
  while in Windows it is generally [PowerShell](https://en.wikipedia.org/wiki/PowerShell **Shell**  
  A shell is the command interpreter that runs inside the terminal. It reads the commands you type, executes them, and displays the results. Common examples include  
  [bash](https://en.wikipedia.org/wiki/Bash_(Unixd  
  [zsh](https://en.wikipedia.org/wiki/MacOS).

You can see the relationship between a human user, a terminal, a shell, and a computer as two people talking in a room.  
One person is you, the other is the computer, and the room where you are talking is the terminal—a common space where you can directly interact.  
The shell is the common set of rules you both understand: your shared language.  
Depending on the operating system and your preferences, you can “speak” with a computer in different rooms and in different languages.  
Generally, *terminal* and *shell* are used interchangeably, although they are not the same.

Before modern graphical interfaces, most interactions with computers happened through a terminal with a keyboard.  
Today most people use a mouse or trackpad to click on windows and files—this is conceptually equivalent to issuing terminal commands.

In scientific contexts, terminal usage is still very common.  
Scientific tools are often built for command‑line use, and supercomputers (such as Baobab at the University of Geneva) provide only terminal access.  
Most importantly, using a terminal—especially in a [Unix](https://en.wikipedia.org/wiki/MacOS)/[Unix‑like](https://en.wikipedia.org/wiki/MacOS) environment like macOS or Ubuntu—provides a universal interface that works the same on laptops, servers, and HPC systems.

This document is a brief introduction to common terminal vocabulary: we will learn **commands** for a given **shell**.  
We focus on commands available in both `bash` and `zsh`, two popular shells.

---

# Part 1 — Getting Comfortable

Start by opening your terminal apps.  
Since this guide focuses on `bash` and `zsh`:

- macOS users → open **Terminal.app**  
- Linux users → open the terminal distributed with your OS  
- Windows users → **do not** use PowerShell (different language)  
  Instead open [Windows Subsystem for Linux (WSL)](https://en.wikipedia.org/wiki/MacOS).

If you set up WSL with us, the first message you see when opening it should indicate you're inside an Ubuntu environment.  
Formally, this is equivalent to being on a Linux system and using its terminal (WSL has no graphical interface).

---

## **Exercise 1: Where Am I?**  
**Goal:** Learn to check the current directory.  
Commands: `pwd`, `ls`

Without a graphical interface, you cannot see where you are located (Desktop? Downloads? etc.).  
You can check using:

a) Run `pwd`  
b) Run `ls`  
c) Write down your working directory and the items inside it

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
