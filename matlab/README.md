# heat1d
Thermal model for planetary applications

## MATLAB MEX-version
This is the MATLAB `MEX`-file implementation of the thermal code. It is based on the original C-code, adapted to be compiled and executed from within the MATLAB environment.

Here are instructions for installing, running, and customizing the model.

First, download the following files and place them in your own directory:

* `heat1d_mex.c` - the core program, including necessary MATLAB MEX-file headers and options, functions, etc.
* `heat1dfun.h` - header file for `heat1dfun.h` containing function prototypes
* `orbitfun.h` - header file for `orbitfun.c`

Open MATLAB. From MATLAB's command line, locate the directory `mydir/` containing these files, and compile the program using:

```
>> mex mydir/heat1d_mex.c
Building with 'Xcode with Clang'.
MEX completed successfully.
```
This should generate a file called `heat1d_mex.mex---`, where `---` is a suffix indicating your operating system.
