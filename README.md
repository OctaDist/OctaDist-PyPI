![Python version][Py-ver-badge]
[![PyPI-Server][PyPI-badge]][PyPI-link]
![Python Wheel][Py-wheel-badge]
![Code size][Code-size]
![Repo size][Repo-size]

[![Travis-CI Test][Travis-badge]][Travis-link]


[Py-ver-badge]: https://img.shields.io/pypi/pyversions/octadist.svg
[PyPI-badge]: https://img.shields.io/pypi/v/octadist.svg
[PyPI-link]: https://pypi.org/project/octadist/
[Py-wheel-badge]: https://img.shields.io/pypi/wheel/octadist.svg
[Code-size]: https://img.shields.io/github/languages/code-size/OctaDist/OctaDist-PyPI.svg
[Repo-size]: https://img.shields.io/github/repo-size/OctaDist/OctaDist-PyPI.svg
[Travis-badge]: https://img.shields.io/travis/OctaDist/OctaDist-PyPI/master.svg
[Travis-link]: https://travis-ci.org/OctaDist/OctaDist-PyPI

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

 

## Standard abilities

OctaDist is computer software for inorganic chemistry and crystallography program.
It is written entirely in Python 3 binding to Tkinter GUI toolkit.
OctaDist can be used for studying the structural distortion in coordination complexes.
With the abilities of OctaDist, you can:

- identify the type of octahedral coordination complexes.
- compute octahedral distortion parameters.
- display 3D molecule and other stuff.
- implement its functions in your or other programs.
- broaden the program with your ideas as needed.

## Documents

User document: [Online manual](https://octadist.github.io/manual.html).

Reference document: [HTML][html] | [PDF][pdf] | [Epub][epub]

[html]: https://octadist-pypi.readthedocs.io/en/latest/ 
[pdf]: https://readthedocs.org/projects/octadist-pypi/downloads/pdf/latest/
[epub]: https://readthedocs.org/projects/octadist/downloads/epub/latest/


## Installing

The easiest way to install OctaDist is `pip` package manager.

```sh
pip install octadist
```

Now you should check if `octadist` package is installed correctly.

```python
import octadist

print(octadist.__version__)     # '2.5.3'
```

## Running the tests

```python
import octadist as oc

# Prepare list of atomic coordinates of octahedral structure:

coord = [[2.298354000, 5.161785000, 7.971898000],   # <- Metal center atom
         [1.885657000, 4.804777000, 6.183726000],   # Ligand atom 1
         [1.747515000, 6.960963000, 7.932784000],   # Ligand atom 2
         [4.094380000, 5.807257000, 7.588689000],   # Ligand atom 3
         [0.539005000, 4.482809000, 8.460004000],   # Ligand atom 4
         [2.812425000, 3.266553000, 8.131637000],   # Ligand atom 5
         [2.886404000, 5.392925000, 9.848966000]]   # Ligand atom 6

zeta = oc.calc_zeta(coord)             # Zeta
delta = oc.calc_delta(coord)           # Delta
sigma = oc.calc_sigma(coord)           # Sigma
theta = oc.calc_theta(coord)           # Theta
```

Example output for computed parameters:

```shell
Computed parameters
-------------------
Zeta  = 0.228072561
Delta = 0.000476251
Sigma = 47.92652837
Theta = 122.6889727
```

Other example scripts and octahedral complexes are available at 
[example-py](https://github.com/OctaDist/OctaDist-PyPI/tree/master/example-py) and 
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

- [Rangsiman Ketkaew](https://sites.google.com/site/rangsiman1993) (Thammasat University, Thailand) <br/>
  - E-mail: rangsiman1993@gmail.com <br/>
- [Yuthana Tantirungrotechai](https://sites.google.com/site/compchem403/people/faculty/yuthana) (Thammasat University, Thailand)
  - E-mail: yt203y@gmail.com
- [David J. Harding](https://www.funtechwu.com/david-j-harding) (Walailak University, Thailand)
  - E-mail: hdavid@mail.wu.ac.th
- [Phimphaka Harding](https://www.funtechwu.com/phimphaka-harding) (Walailak University, Thailand)
  - E-mail: kphimpha@mail.wu.ac.th
- [Mathieu Marchivie](http://www.icmcb-bordeaux.cnrs.fr/spip.php?article562&lang=fr) (Université de Bordeaux, France)
  - E-mail: mathieu.marchivie@icmcb.cnrs.fr
