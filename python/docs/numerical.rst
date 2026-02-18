Numerical Methods
=================

``heat1d`` implements three finite difference schemes for solving the 1-D heat
equation, following the formulation in Hayne et al. (2017), Appendix A2
(Eqs. A13--A20).

Finite Difference Discretization
---------------------------------

The heat equation on a non-uniform grid is discretized using the flux-conservative
form. For interior node :math:`i`, the semi-discrete equation is:

.. math::
   \rho_i c_{p,i} \frac{\partial T_i}{\partial t} =
   p_i K_{i-1} (T_{i-1} - T_i) + q_i K_i (T_{i+1} - T_i)

where the geometric coefficients :math:`p_i` and :math:`q_i` account for the
non-uniform grid spacing:

.. math::
   p_i = \frac{2 \Delta z_i}{\Delta z_{i-1} \Delta z_i (\Delta z_{i-1} + \Delta z_i)}

.. math::
   q_i = \frac{2 \Delta z_{i-1}}{\Delta z_{i-1} \Delta z_i (\Delta z_{i-1} + \Delta z_i)}

Explicit Scheme (Forward Euler)
-------------------------------

The explicit scheme (Eq. A17) updates temperatures using only values from the
current time step:

.. math::
   T_i^{n+1} = T_i^n + \frac{\Delta t}{\rho_i c_{p,i}}
   \left[ \alpha_i T_{i-1}^n - (\alpha_i + \beta_i) T_i^n + \beta_i T_{i+1}^n \right]

where :math:`\alpha_i = p_i K_{i-1}` and :math:`\beta_i = q_i K_i`.

**Stability**: The explicit scheme is conditionally stable, requiring the
CFL condition:

.. math::
   \Delta t \leq F \cdot \min_i \frac{\rho_i c_{p,i} \Delta z_i^2}{K_i}

where :math:`F \leq 0.5` is the Fourier mesh number. For typical Moon parameters,
this limits the time step to :math:`\sim 3000` s (:math:`\sim 830` steps per
lunar day).

Implicit Scheme (Backward Euler)
--------------------------------

The fully implicit scheme evaluates spatial derivatives at the new time level
:math:`n+1`:

.. math::
   -a_i T_{i-1}^{n+1} + (1 + a_i + b_i) T_i^{n+1} - b_i T_{i+1}^{n+1} = T_i^n

where:

.. math::
   a_i = \frac{\Delta t \, p_i K_{i-1}}{\rho_i c_{p,i}}, \quad
   b_i = \frac{\Delta t \, q_i K_i}{\rho_i c_{p,i}}

This forms a tridiagonal system that is solved using the Thomas algorithm
(see below).

**Stability**: The implicit scheme is *unconditionally stable*, meaning the time
step is limited by accuracy rather than stability. This allows much larger time
steps (e.g., :math:`\sim 24` steps per lunar day), providing a :math:`\sim 35\times`
speedup over the explicit scheme.

**Accuracy**: First-order in time (:math:`O(\Delta t)`), same as explicit.

Crank-Nicolson Scheme (Semi-Implicit)
-------------------------------------

The Crank-Nicolson scheme averages the explicit and implicit contributions,
achieving second-order accuracy in time:

.. math::
   -\frac{a_i}{2} T_{i-1}^{n+1} + \left(1 + \frac{a_i + b_i}{2}\right) T_i^{n+1}
   - \frac{b_i}{2} T_{i+1}^{n+1} =
   \frac{a_i}{2} T_{i-1}^n + \left(1 - \frac{a_i + b_i}{2}\right) T_i^n
   + \frac{b_i}{2} T_{i+1}^n

The left-hand side forms a tridiagonal system with half-coefficients, while the
right-hand side uses the old temperature values.

**Stability**: Unconditionally stable, like the fully implicit scheme.

**Accuracy**: Second-order in time (:math:`O(\Delta t^2)`), the most accurate
of the three schemes for a given time step.

**Boundary values**: The Crank-Nicolson scheme requires both old and new
boundary values. The old values (before the BC update) contribute to the explicit
half of the RHS, while the new values contribute to the implicit half of the LHS.

Thomas Algorithm (TDMA)
-----------------------

Both the implicit and Crank-Nicolson schemes produce tridiagonal linear systems
of the form:

.. math::
   \mathbf{A} \mathbf{x} = \mathbf{d}

where :math:`\mathbf{A}` is tridiagonal with sub-diagonal :math:`\ell`, main
diagonal :math:`m`, and super-diagonal :math:`u`.

The Thomas algorithm (Tri-Diagonal Matrix Algorithm) solves this in :math:`O(n)`
time using forward elimination and back-substitution:

**Forward sweep**:

.. math::
   c'_0 = \frac{u_0}{m_0}, \quad
   d'_0 = \frac{d_0}{m_0}

.. math::
   c'_i = \frac{u_i}{m_i - \ell_{i-1} c'_{i-1}}, \quad
   d'_i = \frac{d_i - \ell_{i-1} d'_{i-1}}{m_i - \ell_{i-1} c'_{i-1}}

**Back-substitution**:

.. math::
   x_n = d'_n, \quad
   x_i = d'_i - c'_i x_{i+1}

The algorithm is numerically stable for the heat equation because the coefficient
matrix is always diagonally dominant: :math:`|m_i| = 1 + a_i + b_i > |a_i| + |b_i|`.

Solver Comparison
-----------------

==================  ===========  ===========  ====================
Scheme              Accuracy     Stability    Steps/lunar day
==================  ===========  ===========  ====================
Explicit            O(dt)        Conditional  ~830
Crank-Nicolson      O(dt^2)      Uncondit.    ~24
Implicit            O(dt)        Uncondit.    ~24
==================  ===========  ===========  ====================

For most applications, the Crank-Nicolson scheme is recommended: it offers
second-order accuracy with the same unconditional stability as the fully implicit
scheme.

API Reference
-------------

.. automodule:: heat1d.solvers
   :members:
   :undoc-members:
   :noindex:
