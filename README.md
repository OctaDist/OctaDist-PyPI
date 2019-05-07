![Python version][Py-ver-badge]
[![PyPI-Server][PyPI-badge]][PyPI-link]
[![Conda-Server][Conda-badge]][Conda-link]
![Python Wheel][Py-wheel-badge]
[![Travis-CI Test][Travis-badge]][Travis-link]
![Code size][Code-size]
![Repo size][Repo-size]


[Py-ver-badge]: https://img.shields.io/pypi/pyversions/octadist.svg
[PyPI-badge]: https://img.shields.io/pypi/v/octadist.svg
[PyPI-link]: https://pypi.org/project/octadist/
[Conda-badge]: https://anaconda.org/rangsiman/octadist/badges/version.svg
[Conda-link]: https://anaconda.org/rangsiman/octadist
[Py-wheel-badge]: https://img.shields.io/pypi/wheel/octadist.svg
[Travis-badge]: https://img.shields.io/travis/OctaDist/OctaDist-PyPI/master.svg
[Travis-link]: https://travis-ci.org/OctaDist/OctaDist-PyPI
[Code-size]: https://img.shields.io/github/languages/code-size/OctaDist/OctaDist-PyPI.svg
[Repo-size]: https://img.shields.io/github/repo-size/OctaDist/OctaDist-PyPI.svg

## OctaDist

Octahedral distortion calculator: 
A tool for calculating distortion parameters in coordination complexes. 
https://octadist.github.io/

<p align="center">
   <img alt="molecule" 
   src="https://raw.githubusercontent.com/OctaDist/OctaDist-PyPI/master/images/molecule.png" 
   align=middle 
   width="200pt" />
<p/>

This program is written entirely in Python 3 and tested on PyCharm (Community Edition). 

## Features

OctaDist is designed as a smart tool used for studying the structural distortion in coordinate complexes.
With the abilities of OctaDist, you can:

- identify the type of octahedral coordination complexes.
- compute octahedral distortion parameters.
- display 3D molecule and other stuff.
- implement its functionality in your or other program.


## Getting started

All details of OctaDist is available at [user manual](https://octadist.github.io/manual.html) on the website.


## Installing

- Install the latest version: 
  ```
  pip install octadist
  ```
- Upgrade/downgrade to a specific version: 
  ```
  pip install --upgrade octadist==2.5.1
  ```
  
## Running the test

Prepare two lists of atomic symbols and atomic coordinates, the latter can be stored in array:

```
atom = ['Fe', 'O', 'O', 'N', 'N', 'N', 'N']

coor = [[2.298354000, 5.161785000, 7.971898000],  # <- Metal center atom
        [1.885657000, 4.804777000, 6.183726000],
        [1.747515000, 6.960963000, 7.932784000],
        [4.094380000, 5.807257000, 7.588689000],
        [0.539005000, 4.482809000, 8.460004000],
        [2.812425000, 3.266553000, 8.131637000],
        [2.886404000, 5.392925000, 9.848966000]]
```

Import necessary module for computing the octahedral distortion parameters, called `calc`:

```
from octadist import calc
```

Then calculate all parameters separately, for example:

```
d_bond = octadist.calc_d_bond(coor)         # Bond distance
d_mean = octadist.calc_d_mean(coor)         # Mean distance
zeta = octadist.calc_zeta(coor)             # Zeta
delta = octadist.calc_delta(coor)           # Delta
angle = octadist.calc_bond_angle(coor)      # Bond angle
sigma = octadist.calc_sigma(coor)           # Sigma
theta = octadist.calc_theta(atom, coor)     # Theta
```

or calculate them at once:

```
zeta, delta, sigma, theta = calc.calc_all(coor)
```

Then print all computed parameters:

```
Computed parameters
-------------------
Zeta  = 0.228072561
Delta = 0.000476251
Sigma = 47.92652837
Theta = 122.6889727
```

Example scripts and coordinate files are available at 
[example-py](https://github.com/OctaDist/OctaDist-PyPI/tree/master/example-py) and at
[example-input](https://github.com/OctaDist/OctaDist-PyPI/tree/master/example-input).

## Citation

Please cite this project when you have used OctaDist for scientific publication.

```
OctaDist: A tool for calculating distortion parameters in coordination complexes.
https://octadist.github.io
```

## Bug report

If you found issues in OctaDist, please report us at [here](https://github.com/OctaDist/OctaDist/issues).

## Project team

- [Rangsiman Ketkaew](https://sites.google.com/site/rangsiman1993) (Thammasat University) <br/>
  - E-mail: rangsiman1993@gmail.com <br/>
- [Yuthana Tantirungrotechai](https://sites.google.com/site/compchem403/people/faculty/yuthana) (Thammasat University)
  - E-mail: yt203y@gmail.com
- [David J. Harding](https://www.funtechwu.com/david-j-harding) (Walailak University)
  - E-mail: hdavid@mail.wu.ac.th
- [Phimphaka Harding](https://www.funtechwu.com/phimphaka-harding) (Walailak University, Thailand)
  - E-mail: kphimpha@mail.wu.ac.th
- [Mathieu Marchivie](http://www.icmcb-bordeaux.cnrs.fr/spip.php?article562&lang=fr) (Université de Bordeaux, France)
  - E-mail: mathieu.marchivie@icmcb.cnrs.fr
