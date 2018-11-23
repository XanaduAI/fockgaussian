

#  Fockgaussian

This repository contains the source code used to produce the results presented in *"A faster calculation of  Franck-Condon factors and Fock matrix elements of Gaussian unitaries using loop hafnians"* [arXiv:1811.](https://arxiv.org/abs/1811.).

## Contents

This repository contains:

* A short python library, *fockgaussian* that given the inputs specifying a Gaussian unitary and the photon numbers of the input and output kets calculates the associated probability amplitude using the loop hafnian function of the `hafnian` library.
* A jupyter notebook with examples of the usage of the library and a comparison against the same calculation using the Fock backend of `strawberryfields`

## Requirements

The functions in fockgaussian and the jupyter notebook examples require [strawberryfields](https://github.com/XanaduAI/strawberryfields) and the [hafnian](https://github.com/XanaduAI/hafnian) library.

## Author

Nicolas Quesada

If you are doing any research using this source code, please cite the following paper:

> Nicolas Quesada.  "A faster calculation of  Franck-Condon factors and Fock matrix elements of Gaussian unitaries using loop hafnians" [arXiv:1811.](https://arxiv.org/abs/1811.)

## License

This source code is free and open source, released under the Apache License, Version 2.0.