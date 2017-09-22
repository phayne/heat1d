# heat1d
Thermal model for planetary applications

## MATLAB MEX-version
This is the MATLAB "MEX"-file implementation of the thermal code. It is based on the original C-code, adapted to be compiled and executed from within the MATLAB environment.

Here are instructions for installing, running, and customizing the model.

First, download the following files and place them in your own directory:

* [`lunarThermalModelCustom.m`](https://github.com/phayne/heat1d/blob/master/matlab/lunarThermalModelCustom.m) - wrapper script in MATLAB language
* [`heat1d_mex.c`](https://github.com/phayne/heat1d/blob/master/matlab/heat1d_mex.c) - the core program, including necessary MATLAB MEX-file headers and options, functions, etc.
* [`heat1dfun.h`](https://github.com/phayne/heat1d/blob/master/matlab/heat1dfun.h) - header file for `heat1dfun.h` containing function prototypes
* [`orbitfun.h`](https://github.com/phayne/heat1d/blob/master/matlab/orbitfun.h) - header file for `orbitfun.c`

In order to compile and run the code, you will need a [MATLAB-supported compiler](https://www.mathworks.com/support/compilers.html) installed on your system. There are freely available C-compilers for Windows, Mac OS-X, and Unix. Once the compiler is properly installed, you may need to [change MATLAB's default compiler](https://www.mathworks.com/help/matlab/matlab_external/changing-default-compiler.html) and other `mex` settings.

Next, open MATLAB. From MATLAB's command line, locate the directory `mydir/` containing these files, and compile the program using the following command:

```
>> mex mydir/heat1d_mex.c
Building with 'Xcode with Clang'.
MEX completed successfully.
```
This should generate a file called `heat1d_mex.mex---`, where `---` is a suffix indicating your operating system. The output above may change depending on your operating system and compiler settings.
