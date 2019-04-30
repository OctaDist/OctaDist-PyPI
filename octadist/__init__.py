"""
OctaDist  Copyright (C) 2019  Rangsiman Ketkaew et al.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

----------------------------------------------------------------------

OctaDist version 2.5.1

Octahedral Distortion Analysis
Software website: https://octadist.github.io
Last modified: April 2019

This program was written in Python 3 binding to TkInter GUI platform,
tested on PyCharm (Community Edition) program, and compiled by PyInstaller.

Authors:
Rangsiman Ketkaew            Thammasat University, Thailand    rangsiman1993@gmail.com
Yuthana Tantirungrotechai    Thammasat University, Thailand    yt203y@gmail.com
David J. Harding             Walailak University, Thailand     hdavid@mail.wu.ac.th
Phimphaka Harding            Walailak University, Thailand     kphimpha@mail.wu.ac.th
Mathieu Marchivie            Université de Bordeaux, France    mathieu.marchivie@icmcb.cnrs.fr
"""

program_version = "2.5.1"
program_revision = "April 2019"

name = "octadist"

__all__ = [
    "calc",
    "coord",
    "draw",
    "elements",
    "linear",
    "plane",
    "projection"
    ]
