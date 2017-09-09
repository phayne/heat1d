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

With these files in the same directory, compile the program using a standard compiler such as `gcc`:

`gcc -lm heat1d_moon.c heat1dfun.c orbitfun.c -o heat1d_moon`

The `-lm` flag tells `gcc` to link the mathematics standard library, and the executable program is called `heat1d_moon` in this case. To display the usage, simply run the program with no input arguments:

```
./heat1d_moon

Usage:
  heat1d_moon [lat] [T.I.] [H] [albedo]

    [lat] -- latitude in degrees
    [T.I.] -- thermal inertia at 273 K [SI units] (50 for typical regolith)
    [H] -- H-parameter = scale height of TI increase (0.06 is typical)
    [albedo] -- solar bolometric albedo of surface

```
