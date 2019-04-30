[![Travis-CI Test](https://img.shields.io/travis/OctaDist/OctaDist/master.svg
)](https://travis-ci.org/OctaDist/OctaDist)
[![Code size](https://img.shields.io/github/languages/code-size/OctaDist/OctaDist.svg)](https://github.com/OctaDist/OctaDist)
[![Repo size](https://img.shields.io/github/repo-size/OctaDist/OctaDist.svg)](https://github.com/OctaDist/OctaDist)

## OctaDist
Octahedral distortion calculator: A program for determining the structural distortion of the octahedral complexes. https://octadist.github.io/

<p align="center">
   <img alt="molecule" src="https://raw.githubusercontent.com/OctaDist/OctaDist/master/images/molecule.png" align=middle width="200pt" />
<p/>

This program is written entirely in Python 3 and tested on PyCharm (Community Edition). 

## Installation
- Install the latest version: `pip install octadist`
- Upgrade to the latest version: `pip install --upgrade octadist`
- Upgrade/downgrade to a specific version: `pip install --upgrade octadist==2.5.0`

## Sample usage
Example scripts are available at [here](https://github.com/OctaDist/OctaDist-PyPI/tree/master/example-py).

First of all, you have to import necessary modules for computing the octahedral distortion parameters, called `calc`:

```
from octadist import calc
```

Prepare list (or array) for atomic labels and coordinates:

```
atom = ['Fe', 'O', 'O', 'N', 'N', 'N', 'N']

coor = [[2.298354000, 5.161785000, 7.971898000],
        [1.885657000, 4.804777000, 6.183726000],
        [1.747515000, 6.960963000, 7.932784000],
        [4.094380000, 5.807257000, 7.588689000],
        [0.539005000, 4.482809000, 8.460004000],
        [2.812425000, 3.266553000, 8.131637000],
        [2.886404000, 5.392925000, 9.848966000]]
```

or you can open input file and extract the octahedral structure from input metal complex using a module called `coord`:

```
from octadist import coord, calc
```

For example, input file `full\path\of\your\input\file\Multiple-metals.xyz` 
(other example input files are available at [here](https://github.com/OctaDist/OctaDist-PyPI/tree/master/example-input)):

```
file = r"full\path\of\your\input\file\Multiple-metals.xyz"

atom, coor = coord.extract_octa(file)
```

Then calculate all octahedral parameters

```
d_mean, zeta, delta, sigma, theta = calc.calc_all(atom, coor)
```

and print all computed parameters:

```
All computed parameters
-----------------------
Zeta  = 0.22807256171728651
Delta = 0.0004762517834704151
Sigma = 47.926528379270124
Theta = 122.688972774546
```

## Citation
Please cite this project when you have used OctaDist for scientific publication.

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
