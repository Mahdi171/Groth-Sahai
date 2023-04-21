# Groth-Sahai Proofs

A prototype implementation of Groth and Sahai Proof systems, introduced in "Efficient Non-interactive Proof Systems for Bilinear Groups" paper at EC'08.

## Instruction for Ubuntu 22.04

### Prerequisite Packages

### Install the PBC Stanford library:
Pairbing-Based Cryptography [PBC](https://crypto.stanford.edu/pbc/) library.

```
wget http://crypto.stanford.edu/pbc/files/pbc-0.5.14.tar.gz
tar xf pbc-0.5.14.tar.gz
cd pbc-0.5.14
sudo ./configure.sh
sudo make
sudo make install
sudo make test
```

This [Youtube clip](https://www.youtube.com/watch?v=T0SHn8lMKJA) also gives a detailed instruction on how to set up the PBC library.

### Charm-crypto needs to be installed manually.

- The charm-crypto library should be installed manually from [this repository](https://github.com/JHUISI/charm.git).
Do not use the releases, they do not work. Install from the repo by running the following commands.

```
git clone https://github.com/JHUISI/charm.git
cd charm
sudo ./configure.sh
sudo make
sudo make install
sudo make test
```

Make sure to set the extra `LDFLAGS` so that charm-crypto finds pbc as shown above.

- In case if you are using VS code as the compiler and have installed multiple versions of Python, you might need to change your python interpreter to the compatible version.

```
view --> command palette --> search for Python: Select Interpreter --> choose your compatiable version.
```

- Note that python 3.8 and above seems to be broken for charm-crypto, see [this issue](https://github.com/JHUISI/charm/issues/239).
