# decomotor
Unsupervised Bayesian Ising Approximation for decoding neural-motor dictionaries. Here we present the code corresponding to the paper:
<-to be completed->

The code in this repository follows a method for decoding neural activity in relation to motor behavior, expanding on the approach BIA by Fisher and Mehta (2015).

## Getting Started

### Prerequisites

Tha main code of this project is written in Mathematica v10 (and tested on v11, as well).

The repository also contains a demo in a jupyter notebook (python).

```
Give examples
```

### Installing

The mx files for the functions can work directly. If they do not work, you will need to create the scripts in our local Mathematica, in order to be sure that they will run properly. Go to the folder "src_fun/" and from there run the following command (or open it with Mathematica and run it from there):
```
math -script create_all_functions.m
```

This command creates the functions that are needed by the main programs. Also, there are notebooks for each function showing their contents.

## Running the tests

As a first test, we run the analysis over two dataset, located in the "dataset/" folder. These contain the spiketrains and different behavioral features (such as pitch, amplitude, spectral entropy) for two neurons in the same songbird. The following program will transform these data into binary matrices, and them run the method over them to find relevant words (with and without the behavioral output). Go to the folder "main/", and run (or open it with Mathematica and run it from there):

```
math -script main_data.m
```
This program generates as output files in the "resdump/" folder with details of the codewords found.

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

