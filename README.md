# acoustic_source_localisation repo
AGH eng. thesis project - acoustic source localisation algo overview, tests and implementation for STM32 H563ZI

# matlab_sim

Here I keep MATLAB 2025Rb code of ASL algos simulations:
- Delay-and-Sum  *(in progress)*
- SRP-PHAT
- MUSIC 
- MVDR adaptive beamforming  *(in progress)*
- Maximum Entropy/Maximum Likelihood  *(TODO)*

Used toolboxes:
- Audio Toolbox
- DSP Toolbox
- Phased Array Toolbox
- Communications Toolbox (*to be verified*)

# srp_phat_impl

In this directory you're gonna find my C implementation of SRP-PHAT algorithm's PoC running on PC device with *kissfft* library

How I build it (on Windows):

```
cd {parent_dirs}\srp_phat_impl
mkdir .\build
cd .\build
cmake -G "MinGW Makefiles" -DCMAKE_C_COMPILER=gcc ..
mingw32-make
.\program.exe
```

On Linux you'd probably want to do it the same way, but run `cmake` as usual and use usual `make` command
