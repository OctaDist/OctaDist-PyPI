from octadist import coord, draw

#############
# Example 4 #
#############

# Open and read input file
file = r"C:\Users\Nutt\Desktop\OctaDist-TestFile\Fe\Multiple-metals.xyz"

# Graphical display for octahedral complex
full_atom, full_coor = coord.get_coord(file)
draw.all_atom(full_atom, full_coor)

