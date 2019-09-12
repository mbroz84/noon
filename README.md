# A generator of forward neutrons for ultra-peripheral collisions

### Scope

Repository contains the generator class and two macros showing various usage

### Building / Installation

Generator is a ROOT based class compiled on-the-fly, the macros are ment to run in ROOT. Supported versions are ROOT5 and ROOT6

### Running

root runBreakup.C

This is an afterburner generator an input file with theory curves: one-arm coherent VM photoproduction cross section and invariant mass histogram, are expected. An example file for testing purposes is provided.

See/Use run macro (runBreakup.C) to execute the generator. There are two sub-functions you may use:
- runStandaloneGenerator(): This method by default runs generation of events based on the theory input. 
- computeModelBreakups(): This method using the same theory input as in the previous method will produce the cross section predictions as on Fig.14 in the paper

### Files:

- ENDF_Pb208.txt: Input file with photo-neutron emission spectra for Pb208
- ExampleTheory.root: Example file with theory prediction
- NeutronGenerator.h .cxx: Generator class
- runBreakup.C: ROOT macro to run the generator 
