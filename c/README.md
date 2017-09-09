# heat1d
Thermal model for planetary applications

## C-version
This is the native `C`-language implementation of the thermal code. It is the most computationally efficient version, though perhaps the least user-friendly.

Here are instructions for installing, running, and customizing the model.

First, download the following files and place them in your own directory:

* `heat1d_moon.c` - the core program
* `heat1dfun.c` - subroutines used by the core program
* `heat1dfun.h` - header file for `heat1dfun.h` containing function declarations
* `orbitfun.c` - subroutines for orbital elements and related functions
* `orbitfun.h` - header file for `orbitfun.c`

