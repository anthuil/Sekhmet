
# Release 1, from Version 2.1 (20251022)
# Antoine Thuillier

# This code reads Sekhmet's outputs files and displays graphs of the data.
# The name(s) of the file(s) are arguments of the code.
# The type of graph to print has to be set below ("print_parameter").

# Execution :
# python3 Amon-Ra_Release-1.py Sekhmet_output*
# python3 Amon-Ra_Release-1.py Sekhmet_output_1.csv Sekhmet_output_2.csv

import sys
import numpy as np
import matplotlib.pyplot as plt

print_parameter = 0
# 0 : semi-major axis vs time
# 1 : eccentricity vs time
# 2 : dot(a) vs time

plt.rcParams.update({'font.size': 18}) # Change the base size of all writings in plots.

line_width = 4 # Thickness of the lines in the graphs
Input_file_list = []

for i in range(1, len(sys.argv)):
    Input_file_list.append(sys.argv[i])

# Time [years], M_Star [kg], R_star [m], L_star [W], dot_M_Star [kg/s], Menv_base [m], Renv_base [m], a [m], dot_a [m/s]
time_serie = []
M_Star = []
R_star = []
L_star = []
dot_M_Star = []
Menv_base = []
Renv_base = []
a = [] # semi-major axis [au]
dot_a = []
e = []
Omega_p = [] # Rotation speed of the planet (on itself) [rad/s]

Initial_parameters_planet = []
Initial_parameters_environment = []
Initial_parameters_program = []

Contact_time = [] # Time at which the planets enters into its star.

au = 149597870700   # [m]
R_Sun = 6.957e8     # [m]
M_Sun = 1.9884e30   # [kg]
id_line = 0
iter_file = 0

figure_orbit, ax = plt.subplots()

for Input_file in Input_file_list:
    iter_file += 1

    print("Extracting data from ", Input_file)
    Initial_parameters_planet.append(np.loadtxt(Input_file, dtype=str, max_rows=1))  # will read 1 row
    Initial_parameters_environment.append(np.loadtxt(Input_file, dtype=str, max_rows=1))  # will read 1 row
    Initial_parameters_program.append(np.loadtxt(Input_file, dtype=str, max_rows=1))  # will read 1 row
    datanames = np.loadtxt(Input_file, dtype=str, skiprows=3, max_rows=1) # will read 1 row
    data = np.loadtxt(Input_file, dtype=float, skiprows=5)  # will skip 5 row
    data.shape

    if 0 :
        # Data extraction.
        # 0 : time
        # 1 : M_star
        # 2 : R_star (R_eff)
        # 3 : L_star
        # 4 : dot_M_star
        # 5 : Menv_base
        # 6 : Renv_base
        # 7 : a
        # 8 : dot_a
        # 9 : e
        # 10 : Omega_p
        # 11 : Omega_star
        pass

    time_serie = data[0:, 0]
    R_star = data[0:, 2]
    a = data[0:, 7]
    e = data[0:, 9]
    dot_a = data[0:, 8]

    R_star = R_star * R_Sun / au # Converting to AU.

    # Looking for the survival flag.
    Survival = -1

    file = open(Input_file,'r') # Opening the data file.
    for line in file:
        if line[1:2] == "#": # If the second character of the line is a "#"
            Survival = int(line[14:15]) # Survival parameter | 0 = contact with star | 1 = no contact
            Contact_time.append(float(line[36:54]))
    file.close()

    # Printing the result.
    if print_parameter == 0 : # Semi-major axis
        if Survival == 1:
            ax.plot(time_serie, a, 'g', linewidth=line_width)
        elif Survival == 0:
            ax.plot(time_serie, a, 'r', linewidth=line_width)
        else:
            print("***************************************")
            print("Error in survival assessement.")
            print(Input_file_list[iter_file-1])
            ax.plot(time_serie, a, 'grey', linewidth=line_width)
            print("Survival value : ", Survival)
    elif print_parameter == 1 : # Eccentricity
        ax.plot(time_serie, e, linewidth=line_width)
    elif print_parameter == 2 : # dot(a)
        ax.plot(time_serie, dot_a, linewidth=line_width)
    else :
        print("*** Error in definition of print_parameter ***")

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def au_to_R_Sun(x):
    return x * au / R_Sun

def R_Sun_to_au(x):
    return x * R_Sun / au

if print_parameter == 0 : # Semi-major axis
    ax.fill_between(time_serie, R_star, 0, color='darkred', interpolate=True, label='Star radius') # Remplissage entre 0 et les valeurs y dans data_y1.
    ax.set_xlabel("Time [years]", size=20)
    ax.set_ylabel("Distance [au]", size=20)
    plt.title("Evolution of a planet's orbit while its star ascend the RGB", size=20)

    print("Contact time : ", Contact_time)

    secondary_axis_R_Sun = ax.secondary_yaxis('right', functions=(au_to_R_Sun, R_Sun_to_au))
    secondary_axis_R_Sun.set_ylabel('Distance [R_Sun]', size=20)

elif print_parameter == 1 : # Eccentricity
    ax.set_xlabel("Time [years]", size=20)
    ax.set_ylabel("Eccentricity", size=20)
    plt.title("Evolution of a planet's orbital eccentricity while its star ascend the RGB", size=20)

elif print_parameter == 2 : # dot_a
    ax.set_xlabel("Time [years]", size=20)
    ax.set_ylabel("Variation of separation [m/s²]", size=20)
    #plt.title("Variation of the Evolution of a planet's orbital separation while its star ascend the RGB", size=20)


plt.show()


print("done")





# end
