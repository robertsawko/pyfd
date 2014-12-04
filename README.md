# Python support for CFD

There's a number of calculations I often perform prior to a CFD calculation such as 
 * Colebrook-White friction factor check with single phase pipe
 * Calculation of turbulence intensities for BC specification
 * Wheeler moment inversion algorithm

## General

### Pressure drop

Simple pressure drop calculation based on Colebrook-White. Tested with efunda
calculator.

## Turbulence

### Turbulent BC

## Population balance modelling
### Methods of moments
#### Wheeler's moment inversion

Also known as Modified Chebyshev algorithm (Gautschi 2004). Obtains a set of
weights and abcissas from the set of moments.

#### Realizability check
Checks whether Hadamard matrices have non-negative determinants


