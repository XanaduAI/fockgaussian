

#  Fockgaussian

This repository contains the source code used to produce the results presented in *"
Franck-Condon factors by counting perfect matchings of graphs with loops
"* [J. Chem. Phys. 150, 164113 (2019)](https://aip.scitation.org/doi/10.1063/1.5086387).

## Contents

This repository contains:

* A short python library, *fockgaussian* that given the inputs specifying a Gaussian unitary and the photon numbers of the input and output kets calculates the associated probability amplitude using the loop hafnian function of the [The Walrus](https://github.com/XanaduAI/thewalrus).
* A jupyter notebook with examples of the usage of the library and a comparison against the same calculation using the Fock backend of [Strawberry Fields](https://github.com/XanaduAI/strawberryfields).

## Requirements

The functions in fockgaussian and the jupyter notebook examples require [Strawberry Fields](https://github.com/XanaduAI/strawberryfields) and the [The Walrus](https://github.com/XanaduAI/thewalrus) library.

## Author

Nicolas Quesada

If you are doing any research using this source code, please cite the following paper:

> Nicolas Quesada.  *"
Franck-Condon factors by counting perfect matchings of graphs with loops
"* [J. Chem. Phys. 150, 164113 (2019)](https://aip.scitation.org/doi/10.1063/1.5086387).

## License

This source code is free and open source, released under the Apache License, Version 2.0.
