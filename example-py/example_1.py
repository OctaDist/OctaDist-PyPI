from octadist import calc

#############
# Example 1 #
#############

# Prepare two lists for atomic labels and coordinates (or array for the latter)
atom = ['Fe', 'O', 'O', 'N', 'N', 'N', 'N']

coor = [[2.298354000, 5.161785000, 7.971898000],
        [1.885657000, 4.804777000, 6.183726000],
        [1.747515000, 6.960963000, 7.932784000],
        [4.094380000, 5.807257000, 7.588689000],
        [0.539005000, 4.482809000, 8.460004000],
        [2.812425000, 3.266553000, 8.131637000],
        [2.886404000, 5.392925000, 9.848966000]]

# Calculate all octahedral parameters
d_mean, zeta, delta, sigma, theta = calc.calc_all(atom, coor)

# Show all computed parameters
print("\nEx.1: All computed parameters")
print("-----------------------------")
print("Mean distance =", d_mean)
print("         Zeta =", zeta)
print("        Delta =", delta)
print("        Sigma =", sigma)
print("        Theta =", theta)
