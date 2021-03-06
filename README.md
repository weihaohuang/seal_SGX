## Imported SEAL library in SGX

### SEAL
SEAL [homepage](https://www.microsoft.com/en-us/research/project/simple-encrypted-arithmetic-library/)
> The Simple Encrypted Arithmetic Library (SEAL) is an easy-to-use homomorphic encryption library, developed by researchers in the Cryptography Research Group at Microsoft Research. SEAL is written in C++, and contains .NET wrappers for the public API. It has no external dependencies, so it is easy to compile in many different environments. SEAL uses the Microsoft Research License Agreement, and is free for research use.

The SEAL version is v2.1. Please refer to [v2_withsealv2.3](https://github.com/heartever/seal_SGX/tree/master/v2_withsealv2.3) for SEAL v2.3.
### SGX
SGX [homepage](https://software.intel.com/en-us/sgx)
> This Intel technology is for application developers who are seeking to protect select code and data from disclosure or modification. Intel® SGX makes such protections possible through the use of enclaves, which are protected areas of execution in memory. Application code can be put into an enclave by special instructions and software made available to developers via the Intel SGX SDK. The SDK is a collection of APIs, libraries, documentation, sample source code, and tools that allows software developers to create and debug applications enabled for Intel SGX in C and C++.


### Seal in SGX
perform time consuming operations, such as bootstrapping, inside an SGX enclave.

This project compiles the application (App) and the SGX library (SealEnclaveTest.signed.so). It contains 2 copies of SEAL.  
* App  
Source files: OutEnclaveTest/, Seal_OutEnclave/
* SGX library  
Source files: SealEnclaveTest/, SEAL/
* Configuration file (including encryption parameters)  
Configure.txt

### Usage
Prerequisite:  
* Install Intel SGX linux SDK on a recent platform supporting SGX.  
* A recent GCC/G++ compiler.

One may also test it in simulation mode without SGX support by changing line 5 of ``Makefile`` from ```SGX_MODE ?= HW``` to ```SGX_MODE ?= SIM```. Then:
* make  
* ./App

### Work in progress
Make it work in client-server mode (see the csmod directory mode for more information)
