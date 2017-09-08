# heat1d
Thermal model for planetary science applications

This repository contains source code and examples for the `heat1d` planetary thermal model. In its present form, `heat1d` uses the finite difference method to solve the one-dimensional heat equation:

$$\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}\left(k \frac{\partial T}{\partial z}\right)$$
