# Python support for CFD

[![Join the chat at https://gitter.im/robertsawko/pyfd](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/robertsawko/pyfd?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)  
[![Build Status](https://travis-ci.org/robertsawko/pyfd.svg?branch=master)](https://travis-ci.org/robertsawko/pyfd)

This is a collection of Python modules for routine CFD calculation but also as
a sandbox for more complicated calculations.
 * Colebrook-White friction factor check with single phase pipe
 * Calculation of turbulence intensities for inlet BC specification
 * Wheeler and PD moment inversion algorithms
 * Method of classes method

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


