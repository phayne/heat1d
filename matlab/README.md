# heat1d
Thermal model for planetary applications

## MATLAB MEX-version
This is the MATLAB ["MEX"-file](https://www.mathworks.com/help/matlab/matlab_external/introducing-mex-files.html) implementation of the thermal code. It is based on the original C-code, adapted to be compiled and executed from within the MATLAB environment.

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

To run the model, you can call it directly as follows:

```
>> ks = 7.4e-4; kd = 3.4e-3; rhos = 1100; rhod = 1800; H = 0.068; % Typical thermophysical parameters
>> local_time = 0:0.01:24; % local time with noon = 0
>> albedo = 0.12; % typical surface albedo = 0.12
>> latitude = 0; % equator
>> [t, z] = heat1d_mex(H, rhos, rhod, ks, kd, latitude, albedo, local_time);
```
The main purpose of the wrapper script is to do some error checking, and allow the use of arbitrary z-values. Some examples are provided below.

### Examples

```
>> mex ~/research/thermal_model/code/matlab/heat1d_mex.c
Building with 'Xcode with Clang'.
MEX completed successfully.
>> help lunarThermalModelCustom
  function t = lunarThermalModelCustom( H, rhos, rhod, ks, kd, loctime, ...
                                          depth, latitude, albedo )
  -------------------------------------------------------------------------
  Lunar thermal model for calculating temperatures using standard
  thermophysical properties.
  Inputs [dimensions]:
    H = H-parameter in meters [1]
    rhos = surface density in kg.m-3 [1]
    rhod = deep density in kg.m-3 [1]
    ks = surface conductivity in W.m-1.K-1 [1]
    kd = conductivity at depth in W.m-1.K-1 [1]
    loctime = local time in (lunar) hours PAST NOON [1xn]
    depth = depth in meters [mx1]
    latitude = unsigned latitude in degrees [1]
    albedo = albedo [1]
  Outputs [dimensions]:
    t = temperature in Kelvin [mxn]
 
  Author: Paul O. Hayne (Jet Propulsion Laboratory)
  Date created: September, 2016
  Modified: September, 2017
  
  File dependencies: Requires MATLAB "mex" file containing main thermal
  model, written in C. At the moment, this is called "heat1d_mex.c". There
  are some other dependencies provided in the header of that file.
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
>>
>> % First set all of the thermophysical parameters, boundary conditions, and time:
>> ks = 7.4e-4; kd = 3.4e-3; rhos = 1100; rhod = 1800; H = 0.068; % Typical thermophysical parameters
>> local_time = 0:0.01:24; % local time with noon = 0
>> albedo = 0.12; % typical surface albedo = 0.12
>> z = 0; % surface
>> latitude = 0; % equator
>> % Run the model:
>> t = lunarThermalModelCustom(H, rhos, rhod, ks, kd, local_time, z, latitude, albedo);
>> plot(local_time, t)
>> axis([0, 24, 50, 400])
>> xlabel('Local Time','interpreter','latex')
>> ylabel('Temperature (K)','interpreter','latex')
>> set(gca,'xtick',0:6:24,'xticklabel',{'12:00','18:00','00:00','06:00','12:00'})
```
![Example diurnal temperature curve for the Moon](figures/example_lunar_equatorial_tsurf.png | width=150)
