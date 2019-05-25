# OctaDist  Copyright (C) 2019  Rangsiman Ketkaew et al.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

__author__ = "R. Ketkaew, Y. Tantirungrotechai, D. J. Harding, P. Harding, and M. Marchivie"
__author_full__ = "Rangsiman Ketkaew, Yuthana Tantirungrotechai, David J. Harding, " \
                  "Phimphaka Harding, and Mathieu Marchivie"
__maintainer__ = "Rangsiman Ketkaew"
__copyright__ = "OctaDist  Copyright (C) 2019  Rangsiman Ketkaew et al."
__license__ = "GNU v3"
__credit__ = "OctaDist Development Team"
__email__ = "rangsiman1993@gmail.com"
__version__ = "2.5.3"
__revision__ = "2019.253"
__release__ = "May 2019"
__status__ = "stable"
__title__ = "Octahedral Distortion Calculator"
__description__ = "OctaDist: A tool for calculating distortion parameters in coordination complexes."
__doc__ = "https://octadist.github.io/manual.html"
__website__ = "https://octadist.github.io"
__github__ = "https://github.com/OctaDist/OctaDist"

__all__ = \
    ['calc',
     'coord',
     'draw',
     'elements',
     'linear',
     'plot',
     'projection',
     'tools',
     'calc_d_bond',
     'calc_d_mean',
     'calc_zeta',
     'calc_delta',
     'calc_bond_angle',
     'calc_sigma',
     'calc_theta',
     'count_line',
     'find_metal',
     'extract_file',
     'extract_octa'
     ]

from .src import __src__

# Bring sub-modules in src package to top-level directory
from .src import calc
from .src import coord
from .src import draw
from .src import elements
from .src import linear
from .src import plot
from .src import projection
from .src import tools

# Bring method in sub-modules to top-level directory
from .src.calc import calc_d_bond
from .src.calc import calc_d_mean
from .src.calc import calc_zeta
from .src.calc import calc_delta
from .src.calc import calc_bond_angle
from .src.calc import calc_sigma
from .src.calc import calc_theta
from .src.calc import calc_theta_min
from .src.calc import calc_theta_max

from .src.coord import count_line
from .src.coord import find_metal
from .src.coord import extract_file
from .src.coord import extract_octa

from .src.coord import check_xyz_file
from .src.coord import check_gaussian_file
from .src.coord import check_nwchem_file
from .src.coord import check_orca_file
from .src.coord import check_qchem_file
from .src.coord import get_coord_xyz
from .src.coord import get_coord_gaussian
from .src.coord import get_coord_nwchem
from .src.coord import get_coord_orca
from .src.coord import get_coord_qchem

from .src.draw import all_atom
from .src.draw import all_atom_and_face
from .src.draw import octa
from .src.draw import octa_and_face
from .src.draw import proj_planes
from .src.draw import twisting_faces

from .src.elements import check_atom
from .src.elements import check_radii
from .src.elements import check_color

from .src.linear import norm_vector
from .src.linear import angle_btw_planes
from .src.linear import triangle_area

from .src.plane import find_eq_of_plane

from .src.plot import plot_zeta_sigma
from .src.plot import plot_sigma_theta

from .src.projection import project_atom_onto_line
from .src.projection import project_atom_onto_plane

from .src.tools import find_bonds
from .src.tools import find_faces_octa

from .src.util import calc_fit_plane
from .src.util import plot_fit_plane
from .src.util import calc_rmsd

