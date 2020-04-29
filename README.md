# CryptoSPN
## Implementation of "CryptoSPN: Privacy-preserving Sum-Product Network Inference" (Accepted at [ECAI'20](https://ecai2020.eu/) - pre-print available [here](https://arxiv.org/abs/2002.00801))
#### by *Amos Treiber* ([ENCRYPTO](https://encrypto.de), TU Darmstadt), *Alejandro Molina* ([AI & ML Lab](https://ml-research.github.io/index.html), TU Darmstadt), *Christian Weinert* ([ENCRYPTO](https://encrypto.de), TU Darmstadt), *Thomas Schneider* ([ENCRYPTO](https://encrypto.de), TU Darmstadt), and *Kristian Kersting* ([AI & ML Lab](https://ml-research.github.io/index.html), TU Darmstadt).
----

In this repository, we provide the code for *CryptoSPN*, an extension of [SPFlow](https://github.com/SPFlow/SPFlow) to allow for privacy-preserving SPN inference.

### About
This code provides an interface between [SPFlow](https://github.com/SPFlow/SPFlow) and the [ABY](https://github.com/encryptogroup/ABY) secure computation framework. Details can be found in our [paper](https://arxiv.org/abs/2002.00801). Essentially, it works like [SPFlow](https://github.com/SPFlow/SPFlow)'s C++ functionality, but with an automatic compilation of the exported ABY C++ files into private executables.
This is an experimental research prototype to evaluate benchmarks and *not intended for real-world deployment*. We make no guarantees about its security and correctness.

### Requirements
This code requires all [SPFlow](https://github.com/SPFlow/SPFlow) and [ABY](https://github.com/encryptogroup/ABY) requirements.

### Installation
1. Install [SPFlow](https://github.com/SPFlow/SPFlow) and download *and build* [ABY](https://github.com/encryptogroup/ABY).
2. Enter the paths to your CryptoSPN and [ABY](https://github.com/encryptogroup/ABY) instances in `CryptoSPN/Constants.py`
3. In your CryptoSPN directory, add executable permission to the compile script: `chmod +x compiling/compile.sh`

### Usage
1. Call `spn_to_aby_file(spn)` on a SPFlow `spn` to create an ABY `.cpp` file. CryptoSPN also provides an interface to automatically compile your SPN into an ABY executable by calling `spn_to_aby_exec`. You can specify filenames, precision, and whether to obliviously select client RVs via an oblivious selection network.
2. The created executable can be called like, e.g., `./example -r 0 -a 127.0.0.1 -b 64 -i 50 -f "all_data.txt"`. The following parameters can be set:
```
 -r [Role: 0/1, required]
 -n [Number of parallel operation elements, optional]
 -b [Bit-length, default 32, optional]
 -s [Symmetric Security Bits, default: 128, optional]
 -a [IP-address, default: localhost, optional]
 -p [Port, default: 7766, optional]
 -t [Single test (leave out for all operations), default: off, optional]
 -y [Type of Sharing used, 0: S_BOOL, 1: S_YAO, default: 1, optional]
 -i [Number of iterations for evaluation, default: 1, optional]
 -f [Input file of user RVs, default: all_data.txt, optional]
```
3. See `example.py` for more details.

### Restrictions
1. Supported leave types as of now are: `Bernoulli, Histogram, Gaussian, Poisson`. Note that for Histogram the domain size is leaked.
2. Utility is currently restricted to one bottom-up pass of the network.

### Acknowledgements
This code includes the oblivious selection network ABY implementation of [https://github.com/encryptogroup/PDTE](https://github.com/encryptogroup/PDTE) by Masoud Naderpour.
