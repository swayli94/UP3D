# NLR7301 two-element airfoil

The NLR7301 is a supercritical airfoil which has been modified to a two-element configuration with a non-retractable flap.
This configuration and the measurements were designed for the purpose of CFD validation[1].
The measurements were taken in 1979 and includes detailed pressure distributions, 
transition onset locations and some boundary layer velocity profiles for the angles of attack 6° and 13.1°.
These measurements have been used extensively for another viscid-inviscid interaction method[2]
as well as Navier Stokes codes, e.g.[3]. There is a single measurement set available,
at a Reynolds Number of 2.51E6 and a Mach Number of .185, comprising a set of 16 lift and 3 drag values.

reference: https://github.com/mranneberg/viiflow-examples

(https://github.com/mranneberg/viiflow-examples/blob/master/NLR7301%202-Element%20Airfoil/NLR-7301.ipynb)

The authors note an offset in the change in distance between flap and main airfoil,
as well as a change in twist on the flap. This change is applied below,
by shifting and rotating the flap. Such a modification was also used in the analysis in [6].

```python
# Rotate flap
center = np.r_[0.94,-0.011]
def rotate(deg,center,x):
    c, s = np.cos(deg*np.pi/180), np.sin(deg*np.pi/180)
    A = np.array([[c, s], [-s, c]])
    # The following transposing is done, because numpy subtracts 1x2 arrays from a Nx2 array, 
    # but not 2x1 arrays form an 2xN array.
    return ((A@((x.T-center.T).T)).T+center.T).T

FLAP = rotate(-.25,center,FLAP0)
FLAP[1,:] += 0.0025
```

All calculations have been performed with a Reynolds Number of 2.51E6 and a Mach Number of 0 or 0.185
using a Karman-Tsien correction for the pressure and lift.

[1] B. van den Berg and B. Oskam. Boundary layer measurements on a two-dimensional wing with flap and a comparison with calculations. NLR MP 79034 U

[2] Cebeci, Tuncer, Eric Besnard, and Hsun Chen. Calculation of multielement airfoil flows, including flap wells. 34th Aerospace Sciences Meeting and Exhibit. 1996.

[3] Schwamborn, Dieter, et al. Development of the DLR tau-code for aerospace applications. Proceedings of the International Conference on Aerospace Science and Technology. Bangalore, India: National Aerospace Laboratories, 2008.

[4] Guo, Chuanliang. Effects of turbulence modelling on the analysis and optimisation of high-lift configurations. Master Thesis, Cranfield University

[5] Van Ingen, J. L. The eN method for transition prediction: historical review of work at TU Delft. AIAA, 2008.

[6] Godin, P., D. W. Zingg, and T. E. Nelson. High-lift aerodynamic computations with one-and two-equation turbulence models. AIAA journal 35.2 (1997): 237-243.

[7] Haase, W. et al. ECARP - European Computational Aerodynamics Research Projects: Validation of CFD Codes and Assessment of Turbulence Models. Notes on Numerical Fluid Mechanics, Vol. 58, 1997.
