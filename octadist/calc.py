"""
OctaDist  Copyright (C) 2019  Rangsiman Ketkaew et al.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
"""

from octadist import plane, projection, linear

import numpy as np


def calc_d_mean(c_octa):
    """
    Calculate mean distance parameter (in Angstrom)

    :param c_octa: array - coordinate of octahedral structure

    :return d_mean: float - mean distance
    :return bond_dist: list - individual bond distance
    """
    bond_dist = []
    for i in range(1, 7):
        distance = linear.distance_between(c_octa[0], c_octa[i])
        bond_dist.append(distance)

    i = 0
    sum_distance = 0
    while i < 6:
        sum_distance += bond_dist[i]
        i += 1

    d_mean = sum_distance / 6

    return d_mean, bond_dist


def calc_zeta(d_mean, bond_dist):
    """
    Calculate Zeta parameter (in Angstrom)

         6
    ζ = sum(|dist_i - d_mean|)
        i=1

    Ref: Phys. Rev. B 85, 064114

    :param d_mean: float - mean distance

    :param bond_dist: list - individual bond distance
    :return zeta: float - zeta parameter in degree
    """
    diff_dist = []
    for i in range(6):
        diff_dist.append(abs(bond_dist[i] - d_mean))

    zeta = sum(diff_dist[i] for i in range(6))

    return zeta


def calc_delta(d_mean, bond_dist):
    """
    Calculate Delta parameter
                               2
          1         / d_i - d \
    Δ =  --- * sum | -------- |
          6        \    d    /

    where d_i is individual M-X distance and d is mean M-X distance.
    Ref: DOI: 10.1107/S0108768103026661  Acta Cryst. (2004). B60, 10-20

    :param d_mean: float - mean distance
    :param bond_dist: list - individual bond distance

    :return delta: float - delta parameter (unitless quantity)
    """
    delta = 0
    for i in range(6):
        diff_dist = (bond_dist[i] - d_mean) / d_mean
        delta += (diff_dist * diff_dist) / 6

    return delta


def calc_sigma(c_octa):
    """
    Calculate Sigma parameter

          12
    Σ = sigma < 90 - angle_i >
         i=1

    where angle_i in an unique cis angle
    Ref: J. K. McCusker et al. Inorg. Chem. 1996, 35, 2100.

    :param c_octa: array - coordinate of octahedral structure

    :return sigma: float - sigma parameter in degree
    :return angle_sigma: list - list of 12 unique angles
    """
    angle_sigma = []
    k = 1
    for i in range(1, 7):
        for j in range(i + 1, 7):
            angle = linear.angle_btw_3vec(c_octa[0], c_octa[i], c_octa[j])
            angle_sigma.append(angle)
            k += 1

    # Sort the angle from the lowest to the highest
    angle_sigma.sort()

    # Remove 3 trans angles (last three angles)
    angle_sigma = angle_sigma[:-3]

    max_angle = angle_sigma[-1]

    sigma = 0
    for i in range(len(angle_sigma)):
        sigma = abs(90.0 - angle_sigma[i]) + sigma

    return sigma, max_angle


def calc_theta(a_octa, c_octa, max_angle):
    """
    Calculate Theta parameter

          24
    Θ = sigma < 60 - angle_i >
         i=1

    where angle_i is an unique angle between two vectors of two twisting face.
    Ref: M. Marchivie et al. Acta Crystal-logr. Sect. B Struct. Sci. 2005, 61, 25.

    :param a_octa: list - atom labels of octahedral structure
    :param c_octa: array - coordinate of octahedral structure
    :param max_angle: float - maximum individual cis angle

    :return theta_mean: float - mean Theta value
    :return face_data: list - atomic labels and coordinates of 8 faces
    """
    labels = list(a_octa[1:])
    ligands = list(c_octa[1:])

    _M = c_octa[0]
    N1 = c_octa[1]
    N2 = c_octa[2]
    N3 = c_octa[3]
    N4 = c_octa[4]
    N5 = c_octa[5]
    N6 = c_octa[6]

    # vector from metal to each ligand atom
    _MN1 = N1 - _M
    _MN2 = N2 - _M
    _MN3 = N3 - _M
    _MN4 = N4 - _M
    _MN5 = N5 - _M
    _MN6 = N6 - _M

    ligands_vec = [_MN1, _MN2, _MN3, _MN4, _MN5, _MN6]

    ###########################################
    # Determine the order of atoms in complex #
    ###########################################

    def_change = 6
    for n in range(6):  # This loop is used to identify which N is in line with N1
        test = linear.angle_btw_2vec(ligands_vec[0], ligands_vec[n])
        if test > (max_angle - 1):
            def_change = n

    test_max = 0
    for n in range(6):
        test = linear.angle_btw_2vec(ligands_vec[0], ligands_vec[n])
        if test > test_max:
            test_max = test
            new_change = n

    if def_change != new_change:
        def_change = new_change

    tp = ligands[4]
    ligands[4] = ligands[def_change]
    ligands[def_change] = tp

    tp = labels[4]
    labels[4] = labels[def_change]
    labels[def_change] = tp

    N1 = ligands[0]
    N2 = ligands[1]
    N3 = ligands[2]
    N4 = ligands[3]
    N5 = ligands[4]
    N6 = ligands[5]

    # updated vector from metal to each ligand atom
    _MN1 = N1 - _M
    _MN2 = N2 - _M
    _MN3 = N3 - _M
    _MN4 = N4 - _M
    _MN5 = N5 - _M
    _MN6 = N6 - _M

    # update the ligands_vec list that contains M-N vector
    ligands_vec = [_MN1, _MN2, _MN3, _MN4, _MN5, _MN6]

    for n in range(6):  # This loop is used to identify which N is in line with N2
        test = linear.angle_btw_2vec(ligands_vec[1], ligands_vec[n])
        if test > (max_angle - 1):
            def_change = n

    test_max = 0
    for n in range(6):  # This loop is used to identify which N is in line with N2
        test = linear.angle_btw_2vec(ligands_vec[1], ligands_vec[n])
        if test > test_max:
            test_max = test
            new_change = n

    if def_change != new_change:
        def_change = new_change

    tp = ligands[5]
    ligands[5] = ligands[def_change]  # swapping of the atom (n+1) just identified above with N6
    ligands[def_change] = tp

    tp = labels[5]
    labels[5] = labels[def_change]
    labels[def_change] = tp

    # the new numbering of atoms is stored into the N1 - N6 lists
    N1 = ligands[0]
    N2 = ligands[1]
    N3 = ligands[2]
    N4 = ligands[3]
    N5 = ligands[4]
    N6 = ligands[5]

    # updated vector from metal to each ligand atom
    _MN1 = N1 - _M
    _MN2 = N2 - _M
    _MN3 = N3 - _M
    _MN4 = N4 - _M
    _MN5 = N5 - _M
    _MN6 = N6 - _M

    ligands_vec = [_MN1, _MN2, _MN3, _MN4, _MN5, _MN6]

    for n in range(6):  # This loop is used to identify which N is in line with N3
        test = linear.angle_btw_2vec(ligands_vec[2], ligands_vec[n])
        if test > (max_angle - 1):
            def_change = n

    test_max = 0
    for n in range(6):  # This loop is used to identify which N is in line with N3
        test = linear.angle_btw_2vec(ligands_vec[2], ligands_vec[n])
        if test > test_max:
            test_max = test
            new_change = n

    if def_change != new_change:
        def_change = new_change

    tp = ligands[3]
    ligands[3] = ligands[def_change]  # swapping of the atom (n+1) just identified above with N4
    ligands[def_change] = tp

    tp = labels[3]
    labels[3] = labels[def_change]
    labels[def_change] = tp

    # the new numbering of atoms is stored into the N1 - N6 lists.
    N1 = ligands[0]
    N2 = ligands[1]
    N3 = ligands[2]
    N4 = ligands[3]
    N5 = ligands[4]
    N6 = ligands[5]

    updated_lig = [N1, N2, N3, N4, N5, N6]

    #####################################################
    # Calculate the Theta parameter ans its derivatives #
    #####################################################

    eqOfPlane = []
    indiTheta = []
    allTheta = []  # list gathering the 6 theta angles

    # loop over 8 faces
    for proj in range(8):
        # Find the equation of the plane
        a, b, c, d = plane.find_eq_of_plane(N1, N2, N3)
        eqOfPlane.append([a, b, c, d])

        # Calculation of projection of M, N4, N5 and N6 onto the plane defined by N1, N2, N3
        _MP = projection.project_atom_onto_plane(_M, a, b, c, d)
        N4P = projection.project_atom_onto_plane(N4, a, b, c, d)
        N5P = projection.project_atom_onto_plane(N5, a, b, c, d)
        N6P = projection.project_atom_onto_plane(N6, a, b, c, d)

        # calculate vector from projected metal to projected ligand
        VTh1 = N1 - _MP
        VTh2 = N2 - _MP
        VTh3 = N3 - _MP
        VTh4 = N4P - _MP
        VTh5 = N5P - _MP
        VTh6 = N6P - _MP

        # calculation of 6 theta angles for 1 projection
        a12 = linear.angle_btw_2vec(VTh1, VTh2)
        a13 = linear.angle_btw_2vec(VTh1, VTh3)
        if a12 < a13:
            crossDirect = np.cross(VTh1, VTh2)
        else:
            crossDirect = np.cross(VTh3, VTh1)

        theta1 = linear.angle_btw_2vec_sign(VTh1, VTh4, crossDirect)
        theta2 = linear.angle_btw_2vec_sign(VTh4, VTh2, crossDirect)
        theta3 = linear.angle_btw_2vec_sign(VTh2, VTh5, crossDirect)
        theta4 = linear.angle_btw_2vec_sign(VTh5, VTh3, crossDirect)
        theta5 = linear.angle_btw_2vec_sign(VTh3, VTh6, crossDirect)
        theta6 = linear.angle_btw_2vec_sign(VTh6, VTh1, crossDirect)

        indiTheta.append([theta1, theta2, theta3, theta4, theta5, theta6])

        # sum individual Theta for 1 projection plane
        sumTheta = abs(theta1 - 60) + abs(theta2 - 60) + abs(theta3 - 60) + \
                   abs(theta4 - 60) + abs(theta5 - 60) + abs(theta6 - 60)

        # update Theta into allTheta list
        allTheta.append(sumTheta)

        tp = N2
        N2 = N4
        N4 = N6
        N6 = N3
        N3 = tp

        tp = labels[1]
        labels[1] = labels[3]
        labels[3] = labels[5]
        labels[5] = labels[2]
        labels[2] = tp

        # if the variable proj = 3, Octahedral face permutation face N1N2N3 switches to N1N4N2
        # then to N1N6N4 then N1N3N6 then back to N1N2N3
        if proj == 3:
            tp = N1
            N1 = N5
            N5 = tp
            tp = N2
            N2 = N6
            N6 = tp
            tp = N3
            N3 = N4
            N4 = tp
            tp = labels[0]
            labels[0] = labels[4]
            labels[4] = tp
            tp = labels[1]
            labels[1] = labels[5]
            labels[5] = tp
            tp = labels[2]
            labels[2] = labels[3]
            labels[3] = tp

        # End of the loop that calculate the 8 projections.

    ##################################
    # Summary of computed parameters #
    ##################################

    theta_mean = 0
    for n in range(8):
        theta_mean += allTheta[n]

    theta_mean = theta_mean / 2

    # this list contains the sorted values of Theta from min to max
    sorted_Theta = sorted(allTheta)

    theta_min = 0
    for n in range(4):
        theta_min += sorted_Theta[n]

    theta_max = 0
    for n in range(4, 8):
        theta_max += sorted_Theta[n]

    NewMinTheta = []
    NewMaxTheta = []
    for n in range(4):
        if allTheta[n] < allTheta[n + 4]:
            NewMinTheta.append(allTheta[n])
            NewMaxTheta.append(allTheta[n + 4])
        else:
            NewMaxTheta.append(allTheta[n])
            NewMinTheta.append(allTheta[n + 4])

    min_theta_new = 0
    max_theta_new = 0
    for n in range(4):
        min_theta_new += NewMinTheta[n]
        max_theta_new += NewMaxTheta[n]

    return theta_mean


def calc_all(a_octa, c_octa):
    """
    Calculate D_mean, Zeta, Delta, Sigma, and Theta parameters

    :param a_octa: list - atom labels of octahedral structure
    :param c_octa: array - coordinate of octahedral structure

    :return all_zeta: list - computed zeta
    :return all_delta: list - computed delta
    :return all_sigma: list - computed sigma
    :return all_theta: list - computed theta
    :return all_face: list - 8 faces of octahedral structures
    """
    # Convert tuple or list to array
    c_octa = np.asarray(c_octa)

    # Calculate distortion parameters
    d_mean, bond_dist = calc_d_mean(c_octa)
    zeta = calc_zeta(d_mean, bond_dist)
    delta = calc_delta(d_mean, bond_dist)
    sigma, max_angle = calc_sigma(c_octa)
    theta = calc_theta(a_octa, c_octa, max_angle)

    return d_mean, zeta, delta, sigma, theta

