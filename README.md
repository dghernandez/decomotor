# uBIA
Unsupervised Bayesian Ising Approximation (uBIA) for decoding neural-motor dictionaries. Here we present the code corresponding to the paper:
<-to be completed->

The code in this repository follows a method for decoding neural activity in relation to motor behavior, expanding on the approach BIA by Fisher and Mehta (2015).

## Getting Started

### Prerequisites

Tha main code of this project is written in Mathematica v10 (and tested on v11, as well).

The repository also contains a demo in python with examples in jupyter notebooks (see python folder).

### Loading the Package

Place the file "asPackage/uBIAmotor.m" in your Package folder (see $Path variable). Then from any notebook where you need the package, just type:
```
<<uBIAmotor`
```
or 
```
Needs["uBIAmotor`"]
```

## Decoding a dictionary

First, we run the analysis over two dataset, located in the "dataset/" folder. These contain the spiketrains and different behavioral features (such as pitch, amplitude, spectral entropy) for two neurons in the same songbird. The following program will transform these data into binary matrices, and then it will run the method over them to find relevant words (with and without the behavioral output). Go to the folder "main/", and run (or open it with Mathematica and run it from there*):

```
math -script uBIA_main_data.m
```
This program generates as output several files in the "resdump/" folder with details of the codewords found.

We have also included a program for testing the method on synthetic data (as described in the appendix of the paper).
Go to the folder "main/", and run (or open it with Mathematica and run it from there*):

```
math -script uBIA_main_synthetic.m
```
This program generates as output the performance of the method and append the result to a file in the "resdump/res_synth/" folder. Within this program, there are options to generate samples from different distributions.

*If you are running these programs from Mathematica directly, make sure to edit the first lines in these files, such that you are working in the proper directory.

## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details

