#

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import sys
from sympy import *
import sympy as sym
import os
from itertools import chain

# We set them to be global:
#global E0_C_I,\
#       V0_C_I,\
#       B0_C_I,\
#       B0_prime_C_I,\
#       E0_14,\
#       V0_14,\
#       B0_14,\
#       B0_prime_14 
#

# Intial candidates for fit, per FU: - thus, the E vs V input data has to be per FU
E0_init = -941.510817926696   
V0_init = 63.54960592453  
B0_init = 76.3746233515232  
B0_prime_init = 4.05340727164527 

def BM(V, E0, V0, B0, B0_prime):
        return  E0+ (2.293710449E+17)*(1E-21)*( (9.0/16.0)*(V0*B0) * (  (((V0/V)**(2.0/3.0)-1.0)**3.0)*B0_prime  + ((V0/V)**(2.0/3.0)-1)**2  *  (6.0-4.0*(V0/V)**(2.0/3.0))  ))


def P(V, E0, V0, B0, B0_prime):
    f0=(3.0/2.0)*B0
    f1=((V0/V)**(7.0/3.0))-((V0/V)**(5.0/3.0))
    f2=((V0/V)**(2.0/3.0))-1
    pressure= f0*f1*(1+(3.0/4.0)*(B0_prime-4)*f2)
    return pressure

def H(V, E0, V0, B0, B0_prime):
     return BM(V, E0, V0, B0, B0_prime) + P(V, E0, V0, B0, B0_prime) * V

filefolder_for_E_and_P_vs_V = '/home/david/Trabajo/structures/Calcite_I_and_II/Scaling_more_volumes/scaling_from_117.743646/B3LYP__SHRINK_8_8__bipolar_18_18__TOLINTEG_8_18__XXLGRID_TOLDEE_8'

# Calcite I (Red triangles):
V_not_p_f_unit_C_I, E_not_p_f_unit_C_I = np.loadtxt(os.path.join(filefolder_for_E_and_P_vs_V, 'calcite_I.dat'), skiprows = 1).T

# Calcite II - Trapped (Green triangles):
V_not_p_f_unit_C_II, E_not_p_f_unit_C_II = np.loadtxt(os.path.join(filefolder_for_E_and_P_vs_V, 'calcite_II_trapped.dat'), skiprows = 1).T

# 14 (Empty grey triangles):
V_14_not_p_f_unit, E_14_not_p_f_unit = np.loadtxt(os.path.join(filefolder_for_E_and_P_vs_V, 'E_vs_V_Calcite_II_including_additional_vols.dat'), skiprows = 1).T
print 'V_14_not_p_f_unit = ', V_14_not_p_f_unit
print 'len(V_14_not_p_f_unit) = ', len(V_14_not_p_f_unit)
# In order to avoid last three
#V_14_not_p_f_unit = V_14_not_p_f_unit[:-5]
#print 'V_14_not_p_f_unit = ', V_14_not_p_f_unit
#print 'len(V_14_not_p_f_unit) = ',len(V_14_not_p_f_unit)
#E_14_not_p_f_unit = E_14_not_p_f_unit[:-5]
#print 'E_14_not_p_f_unit = ', E_14_not_p_f_unit
#sys.exit()

# Calcite II - Trapped 0.98 (Green triangles):
V_not_p_f_unit_C_II_0_98, E_not_p_f_unit_C_II_0_98 = np.loadtxt(os.path.join(filefolder_for_E_and_P_vs_V, 'calcite_II_trapped_0_98.dat'), skiprows = 1).T

# Calcite II - Trapped 0.87 (Green triangles):
V_not_p_f_unit_C_II_0_87, E_not_p_f_unit_C_II_0_87 = np.loadtxt(os.path.join(filefolder_for_E_and_P_vs_V, 'calcite_II_trapped_0_87.dat'), skiprows = 1).T

# If the data is not per f unit, do this:
nFU_C_I = 2.0
nFU_C_II = 4.0
E_C_I = E_not_p_f_unit_C_I/nFU_C_I
V_C_I = V_not_p_f_unit_C_I/nFU_C_I

E_C_II = E_not_p_f_unit_C_II/nFU_C_II
V_C_II = V_not_p_f_unit_C_II/nFU_C_II

E_14 = E_14_not_p_f_unit/nFU_C_II
V_14 = V_14_not_p_f_unit/nFU_C_II

E_C_II_0_98 = E_not_p_f_unit_C_II_0_98/nFU_C_II
V_C_II_0_98 = V_not_p_f_unit_C_II_0_98/nFU_C_II

E_C_II_0_87 = E_not_p_f_unit_C_II_0_87/nFU_C_II
V_C_II_0_87 = V_not_p_f_unit_C_II_0_87/nFU_C_II


######

E0_init_14 =  -941.550233907
V0_init_14 =  72.083297156
B0_init_14 =  16.0099525503
B0_prime_init_14 =  6.15613681931


init_vals = [E0_init, V0_init, B0_init, B0_prime_init]
init_vals_14 = [E0_init_14, V0_init_14, B0_init_14, B0_prime_init_14]

popt_C_I, pcov_C_I = curve_fit(BM, V_C_I, E_C_I, p0=init_vals)
popt_C_II, pcov_C_II = curve_fit(BM, V_C_II, E_C_II, p0=init_vals)
popt_14, pcov_14 = curve_fit(BM, V_14, E_14, p0=init_vals_14)

# !!!!!!!!! BEGIN OF THE SYMBOLIC COMMON TANGENT
## Assigning the parameters:
#E0_C_I = popt_C_I[0]
#V0_C_I = popt_C_I[1]
#B0_C_I = popt_C_I[2]
#B0_prime_C_I = popt_C_I[3]
#
#E0_14 = popt_14[0]
#V0_14 = popt_14[1]
#B0_14 = popt_14[2]
#B0_prime_14 = popt_14[3]
#
#


# Linspace for plotting the fitting curves:
V_C_I_lin = np.linspace(V_C_I[0], V_C_I[-1], 10000)
V_C_II_lin = np.linspace(V_C_II[0], V_C_II[-1], 10000)
V_14_lin = np.linspace(V_14[0], V_14[-1], 10000)

# Plotting the fitting curves:
p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='grey', label='BM fit Calcite I' )
p4, = plt.plot(V_C_II_lin, BM(V_C_II_lin, *popt_C_II), 'k--', label='BM fit Calcite II ("trapped")')
p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='BM fit Calcite II')

# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
p3 = plt.scatter(V_C_II, E_C_II, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)
p5 = plt.scatter(V_14, E_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

fontP = FontProperties()
fontP.set_size('small')

plt.legend((p1, p2, p3, p4, p5, p6), ("Calcite I", "BM fit Calcite I", 'Calcite II ("trapped")', 'BM fit Calcite II ("trapped")', "Calcite II", 'BM fit Calcite II'), prop=fontP)

pressures_per_F_unit_C_I = P(V_C_I, *popt_C_I)
print 'CI = ', popt_C_I
output_array_2 = np.vstack((E_not_p_f_unit_C_I, V_not_p_f_unit_C_I, E_C_I, V_C_I, pressures_per_F_unit_C_I)).T
np.savetxt('Volumes_and_pressures.dat', output_array_2, header="Energy (a.u.) \t Volume (A^3) \t Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressures (GPa)", fmt="%0.13f")

pressures_per_F_unit_14 = P(V_14, *popt_14)
print 'C14 = ', popt_14

plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
plt.ticklabel_format(useOffset=False)

plt.savefig('calcite_I_and_II_all_2_summary_better_plot.pdf', bbox_inches='tight')

# Plotting only CI and 14:
plt.figure()

# Plotting the fitting curves:
p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='black', label='BM fit Calcite I' )
p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='BM fit Calcite II')

# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(V_14, E_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

fontP = FontProperties()
fontP.set_size('small')

plt.legend((p1, p2, p5, p6), ("Calcite I", "BM fit Calcite I", "Calcite II", 'BM fit Calcite II'), prop=fontP)

pressures_per_F_unit_C_I = P(V_C_I, *popt_C_I)
print 'CI = ', popt_C_I
output_array_2 = np.vstack((E_not_p_f_unit_C_I, V_not_p_f_unit_C_I, E_C_I, V_C_I, pressures_per_F_unit_C_I)).T
np.savetxt('Volumes_and_pressures.dat', output_array_2, header="Energy (a.u.) \t Volume (A^3) \t Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressures (GPa)", fmt="%0.13f")

pressures_per_F_unit_14 = P(V_14, *popt_14)
print 'C14 = ', popt_14

plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
plt.ticklabel_format(useOffset=False)

plt.savefig('calcite_I_and_II_all_2_summary_better_plot_only_CI_and_14.pdf', bbox_inches='tight')

plt.figure()
# Plotting the fitting curves:
p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='black', label='BM fit Calcite I' )

# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)

fontP = FontProperties()
fontP.set_size('small')

plt.legend((p1, p2), ("Calcite I", "BM fit Calcite I"), prop=fontP)

plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
plt.ticklabel_format(useOffset=False)
plt.savefig('calcite_I_summary_better_plot.pdf', bbox_inches='tight')
plt.clf()
plt.cla()
plt.close()
plt.figure()

# Plotting the fitting curves:
p4, = plt.plot(V_C_II_lin, BM(V_C_II_lin, *popt_C_II), 'k--', label='BM fit Calcite II ("trapped")')

# Plotting the scattered points: 
p3 = plt.scatter(V_C_II, E_C_II, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)

fontP = FontProperties()
fontP.set_size('small')

plt.legend((p3, p4), ('Calcite II ("trapped")', 'BM fit Calcite II ("trapped")'), prop=fontP)

plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
plt.ticklabel_format(useOffset=False)
plt.savefig('calcite_II_summary_better_plot.pdf', bbox_inches='tight')
plt.clf()
plt.cla()

plt.figure()

# Plotting the scattered points: 
plt.scatter(V_C_II_0_87, E_C_II_0_87, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)

fontP = FontProperties()
fontP.set_size('small')

plt.legend() #(p3, p4), ('Calcite II ("trapped")', 'BM fit Calcite II ("trapped")'), prop=fontP)

plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 0.98)$V_{eq}$", fontsize=10)
plt.ticklabel_format(useOffset=False)
plt.ylim(-941.549, -941.535)
plt.savefig('calcite_II_0_87_summary_better_plot.pdf', bbox_inches='tight')
plt.clf()
plt.cla()


plt.figure()

# Plotting the scattered points: 
plt.scatter(V_C_II_0_98, E_C_II_0_98, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)

fontP = FontProperties()
fontP.set_size('small')

plt.legend() #p3, 'Calcite II ("trapped")' ) #, prop=fontP)
plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$", fontsize=10)
plt.ylim(-941.55, -941.542) # adjust
plt.ticklabel_format(useOffset=False)
plt.savefig('calcite_II_0_98_summary_better_plot.pdf', bbox_inches='tight')
plt.clf()
plt.cla()

# H = E + PV

H_C_I = E_C_I + pressures_per_F_unit_C_I * V_C_I
print 'H_C_I =  ', H_C_I
print 'E_C_I =  ', E_C_I
print 'V_C_I =  ', V_C_I

H_14 = E_14 + pressures_per_F_unit_14 * V_14

output_array_3 = np.vstack((E_C_I, V_C_I, pressures_per_F_unit_C_I, H_C_I)).T
np.savetxt('E_V_P_H__C_I.dat', output_array_3, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressure / F.U. (GPa) \t Enthalpy (a.u.)", fmt="%0.13f") 

output_array_4 = np.vstack((E_14, V_14, pressures_per_F_unit_14, H_14)).T
np.savetxt('E_V_P_H__14.dat', output_array_4, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressure / F.U. (GPa) \t Enthalpy (a.u.)", fmt="%0.13f") 



# Fitting Delta_H;

# Quadratic fit of T=T(P):
#c, d, f = np.polyfit(P1, T1, 2)
fitting_C_I = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 1)
fit_C_I = np.poly1d(fitting_C_I)

fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 1)
fit_14 = np.poly1d(fitting_14)


print """
HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
fit_C_I = """, fit_C_I
print 'fit_14 = ', fit_14

fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 2)
fit = np.poly1d(fitting)

print """
HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
fit = """, fit

fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 3)
fit = np.poly1d(fitting)

print """
HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
fit = """, fit


# If we want the Regression coefficcient:
# Polynomial Regression
def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results

All_in_one_1st_degree = polyfit(pressures_per_F_unit_C_I, H_C_I, 1)
All_in_one_2nd_degree = polyfit(pressures_per_F_unit_C_I, H_C_I, 2)
All_in_one_3rd_degree = polyfit(pressures_per_F_unit_C_I, H_C_I, 3)
All_in_one_4th_degree = polyfit(pressures_per_F_unit_C_I, H_C_I, 4)

print 'All_in_one_1st_degree = ',All_in_one_1st_degree
print 'All_in_one_2nd_degree = ',All_in_one_2nd_degree
print 'All_in_one_3rd_degree = ',All_in_one_3rd_degree
print 'All_in_one_4th_degree = ',All_in_one_4th_degree


# Plotting Delta_H:
#********* 1st degree:
fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 1)
fit = np.poly1d(fitting)

fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 1)
fit_14 = np.poly1d(fitting_14)

fig = plt.figure()

EnergyCI, VolumeCI, PressureCI, EnthalpyCI  = np.loadtxt('./E_V_P_H__C_I.dat', skiprows = 1).T

print 'PressureCI[0] = ', PressureCI[0]
print 'PressureCI[-1] = ', PressureCI[-1]

Energy14, Volume14, Pressure14, Enthalpy14  = np.loadtxt('./E_V_P_H__14.dat', skiprows = 1).T

print 'Pressure14[0] = ', Pressure14[0]
print 'Pressure14[-1] = ', Pressure14[-1]

xp_C_I = np.linspace(PressureCI[-1], PressureCI[0], 100)
xp_14 = np.linspace(Pressure14[-1], Pressure14[0], 100)


# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

p11, = plt.plot(xp_C_I, fit(xp_C_I), "black", label='linear fit')
p55, = plt.plot(xp_14, fit_14(xp_14), "grey", label='linear fit')


fontP = FontProperties()
fontP.set_size('small')

plt.legend((p1, p11, p5, p55), ("Calcite I", 'linear fit', 'Calcite II', 'linear fit'), prop=fontP)

plt.xlabel('P / Formula unit (GPa)')
plt.ylabel('H / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
plt.ticklabel_format(useOffset=False)

print 'fitting = ', fitting
print 'fitting_14 = ', fitting_14
crossing_x = np.roots(fitting - fitting_14)
crossing_y = fit(crossing_x)
print 'crossing_x = ', crossing_x
print 'crossing_y = ', crossing_y
#ax = fig.add_subplot(111)
#ax.annotate('Intersection\nP=6.76689272 GPa\nH = -551.95344408 a.u.', xy=(6.76689272, -551.95344408), xytext=(6.76689272+2.7767, -551.95344408-162.27),
#           arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
#           )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_1st_degree.pdf', bbox_inches='tight')


#********* 2nd degree:
fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 2)
fit = np.poly1d(fitting)

fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 2)
fit_14 = np.poly1d(fitting_14)

fig = plt.figure()
xp_C_I = np.linspace(PressureCI[-1], PressureCI[0], 100)
xp_14 = np.linspace(Pressure14[-1], Pressure14[0], 100)

# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

p11, = plt.plot(xp_C_I, fit(xp_C_I), "black", label='2nd degree pol. fit')
p55, = plt.plot(xp_14, fit_14(xp_14), "grey", label='2nd degree pol. fit')

fontP = FontProperties()
fontP.set_size('small')

plt.legend((p1, p11, p5, p55), ("Calcite I", '2nd degree pol. fit', 'Calcite II', '2nd degree pol. fit'), prop=fontP)

plt.xlabel('P / Formula unit (GPa)')
plt.ylabel('H / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)

plt.ticklabel_format(useOffset=False)
crossing_x = np.roots(fitting - fitting_14)
crossing_y = fit(crossing_x)
print 'crossing_x = ', crossing_x
print 'crossing_y = ', crossing_y
#ax = fig.add_subplot(111)
#ax.annotate('Intersection\nP= 2.60943134 GPa\nH = -779.63314906 a.u.', xy=(2.60943134, -779.63314906), xytext=(2.60943134+2.7767, -779.63314906-162.27),
#           arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
#           )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_2nd_degree.pdf', bbox_inches='tight')


#********* 3rd degree:
fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 3)
fit = np.poly1d(fitting)

fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 3)
fit_14 = np.poly1d(fitting_14)

fig = plt.figure()
xp_C_I = np.linspace(PressureCI[-1], PressureCI[0], 100)
xp_14 = np.linspace(Pressure14[-1], Pressure14[0], 100)


# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

p11, = plt.plot(xp_C_I, fit(xp_C_I), "black", label='3rd degree pol. fit')
p55, = plt.plot(xp_14, fit_14(xp_14), "grey", label='3rd degree pol. fit')

fontP = FontProperties()
fontP.set_size('small')

plt.legend((p1, p11, p5, p55), ("Calcite I", '3rd degree pol. fit', 'Calcite II', '3rd degree pol. fit'), prop=fontP)

plt.xlabel('P / Formula unit (GPa)')
plt.ylabel('H / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)

plt.ticklabel_format(useOffset=False)
crossing_x = np.roots(fitting - fitting_14)
crossing_y = fit(crossing_x)
print 'crossing_x = ', crossing_x
print 'crossing_y = ', crossing_y
#ax = fig.add_subplot(111)
#ax.annotate('Intersection\nP =  [ -5.18698424+43.15109474j\n       -5.18698424-43.15109474j\n      1.50582599 +0.j  ]\nH =  [  811.00634615+1866.26303147j\n      811.00634615-1866.26303147j\n       -846.40692208   +0.j ]', xy=(1.50582599, -846.40692208), xytext=(1.50582599+2.7767, -846.40692208-320.27),
#           arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
#           )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_3rd_degree.pdf', bbox_inches='tight')

#********* 4th degree:
fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 4)
fit = np.poly1d(fitting)

fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 4)
fit_14 = np.poly1d(fitting_14)

fig = plt.figure()
xp_C_I = np.linspace(PressureCI[-1], PressureCI[0], 100)
xp_14 = np.linspace(Pressure14[-1], Pressure14[0], 100)


# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

p11, = plt.plot(xp_C_I, fit(xp_C_I), "black", label='4th degree pol. fit')
p55, = plt.plot(xp_14, fit_14(xp_14), "grey", label='4th degree pol. fit')

fontP = FontProperties()
fontP.set_size('small')

plt.legend((p1, p11, p5, p55), ("Calcite I", '4th degree pol. fit', 'Calcite II', '4th degree pol. fit'), prop=fontP)

plt.xlabel('P / Formula unit (GPa)')
plt.ylabel('H / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)

plt.ticklabel_format(useOffset=False)
crossing_x = np.roots(fitting - fitting_14)
crossing_y = fit(crossing_x)
print 'crossing_x = ', crossing_x
print 'crossing_y = ', crossing_y
#ax = fig.add_subplot(111)
#ax.annotate('Intersection\nP=1.05995113 GPa\nH = -874.29224237 a.u.', xy=(1.05995113, -874.29224237), xytext=(1.05995113+2.7767, -874.29224237-162.27),
#           arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
#           )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_4th_degree.pdf', bbox_inches='tight')


#sys.exit()
# Plotting P vs V ADDITIONAL VOLUMES: E_vs_V_Calcite_II_including_additional_vols.dat
fig = plt.figure()
EnergyCI, VolumeCI, PressureCI, EnthalpyCI  = np.loadtxt('./E_V_P_H__C_I.dat', skiprows = 1).T

print 'PressureCI[0] = ', PressureCI[0]
print 'PressureCI[-1] = ', PressureCI[-1]

Energy14, Volume14, Pressure14, Enthalpy14  = np.loadtxt('./E_V_P_H__14.dat', skiprows = 1).T

print 'Pressure14[0] = ', Pressure14[0]
print 'Pressure14[-1] = ', Pressure14[-1]

xp_C_I = np.linspace(PressureCI[-1], PressureCI[0], 100)
xp_14 = np.linspace(Pressure14[-1], Pressure14[0], 100)


# Plotting the fitting curves:
# Initial (this works):
p2, = plt.plot(V_C_I_lin, P(V_C_I_lin, *popt_C_I), color='black', label='BM fit Calcite I' )

# The fix that worked in H
#p2, = plt.plot(V_C_I, P(V_C_I, *popt_C_I), color='grey', label='H fit Calcite II')

# Initial (this works):
p6, = plt.plot(V_14_lin, P(V_14_lin, *popt_14), 'b', label='BM fit Calcite II')

# The fix that worked in H
#p6, = plt.plot(V_14, P(V_14, *popt_14), 'b', label='BM fit Calcite II')

# Plotting the scattered points: 
p1 = plt.scatter(VolumeCI, PressureCI, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(Volume14, Pressure14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

fontP = FontProperties()
fontP.set_size('small')

plt.legend((p1, p2, p5, p6), ("Calcite I", "BM fit Calcite I", 'Calcite II', 'BM fit Calcite II'), prop=fontP)

plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('P (GPa)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
plt.ticklabel_format(useOffset=False)

plt.savefig('calcite_I_and_II_all_2_summary_better_plot_P_vs_V.pdf', bbox_inches='tight')

#plt.show()

# Saving into variables:
P_lin_C_I = P(V_C_I_lin, *popt_C_I)
H_lin_C_I = H(V_C_I_lin, *popt_C_I)

P_lin_14 = P(V_14_lin, *popt_14)
H_lin_14 = H(V_14_lin, *popt_14)

print ' P_lin_C_I = ', P_lin_C_I 
print ' type(P_lin_C_I) = ', type(P_lin_C_I) 
print ' H_lin_C_I = ', H_lin_C_I
print ' P_lin_14  = ', P_lin_14 
print ' H_lin_14  = ', H_lin_14 

output_array_1 = np.vstack((P_lin_C_I, H_lin_C_I)).T
np.savetxt('P_lin_C_I__H_lin_C_I.dat', output_array_1, header="P(GPa) \t   H per F unit (a.u)", fmt="%0.13f")

output_array_2 = np.vstack((P_lin_14, H_lin_14)).T
np.savetxt('P_lin_14__H_lin_14.dat', output_array_2, header="P(GPa) \t    H per F unit (a.u)", fmt="%0.13f")



# Plotting and fitting to the exact expression of Delta_H, for checking the data
#********* Exact expression of Delta_H, for checking the data:

fig = plt.figure()

P_lin_C_I, H_lin_C_I  = np.loadtxt('P_lin_C_I__H_lin_C_I.dat', skiprows = 1).T
P_lin_14, H_lin_14 = np.loadtxt('P_lin_14__H_lin_14.dat', skiprows = 1).T

# Plotting the fitting curves:

# IOBE suggestion:
#p2, = plt.plot(P(V_C_I_lin, *popt_C_I), H(V_C_I_lin, *popt_C_I), color='grey', label='H fit Data' )

p2, = plt.plot(P_lin_C_I, H_lin_C_I, color='magenta', label='H fit from Data' )
p6, = plt.plot(P_lin_14, H_lin_14, color='yellow', label='H fit from Data' )

# IOBE suggestion:
#p6, = plt.plot(P(V_14_lin, *popt_14), H(V_14_lin, *popt_14), color='blue', label='H fit Data' )

print 'V_C_I_lin = ', V_C_I_lin
print ' xp_C_I = ', xp_C_I
print ' H(V_C_I_lin, *popt_C_I) = ', H(V_C_I_lin, *popt_C_I)

# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

fontP = FontProperties()
fontP.set_size('small')

plt.legend((p1, p2, p5, p6), ("Calcite I", 'H fit Calcite I', 'Calcite II', 'H fit Calcite II'), prop=fontP)

plt.xlabel('P / Formula unit (GPa)')
plt.ylabel('H / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
plt.ticklabel_format(useOffset=False)


plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_exact_expression_from_data.pdf', bbox_inches='tight')

#plt.show()
#sys.exit()



#P_lin_C_I, H_lin_C_I = np.loadtxt('./P_lin_C_I__H_lin_C_I.dat').T
#P_lin_14, H_lin_14 = np.loadtxt('./P_lin_14__H_lin_14.dat').T

print 'Performing the collisions program....'

def within_tolerance(p1, p2):
    # 1000 points:
#   tol = 2e-4       # Tests: 1e-4 (no collisions), 1e-3 (two collisions): collisions in calcite 1 =  [(5.2810550218375001, -622.97533180077869), (2.4948809013169, -786.19369268484559)] #collisions in calcite 2 =  [(5.3259729117925003, -622.9745825047172), (2.4889672325627998, -786.1933317892765)]

    # 10000 points:
#   tol = 1e-5    # 1e-6 no collisions. 1e-5 (two collisions): collisions in calcite 1 =  [(7.7584638898975999, -485.33130707345208), (6.8467697399844996, -535.24452210341201)] #collisions in calcite 2 =  [(8.1461924560147008, -485.33130018244401), (7.0184016533328002, -535.24451499780525)]

    # 10000 points and imposing both tolerances:
    tol = 1e-3   #1e-4: no collisions   #1e-3 (two collisions): collisions in calcite 1 =  [(3.5382670841772002, -723.91821227064997), (3.5287623836802, -724.47900542569573), (3.5192624232179002, -725.03963525831909), (3.4491097907470998, -729.18322571800911)] #collisions in calcite 2 =  [(3.5385176471450999, -723.91895811704467), (3.5289151962848, -724.47931966566648), (3.5193081152676999, -725.04013049379648), (3.4484235192196002, -729.18352013774256)]

 
    P_lin_C_I, H_lin_C_I = p1
    P_lin_14, H_lin_14 = p2

    return abs(H_lin_C_I - H_lin_14) < tol and abs(P_lin_C_I - P_lin_14) < tol


points_1 = list(zip(P_lin_C_I, H_lin_C_I))
points_2 = list(zip(P_lin_14, H_lin_14))


collisions_2 = []

for p1 in points_1:
    matches = [p2 for p2 in points_2 if within_tolerance(p1, p2)]
    collisions_2.append(matches)

collisions_1 = []

for p2 in points_2:
    matches = [p1 for p1 in points_1 if within_tolerance(p1, p2)]
    collisions_1.append(matches)

print 'collisions_1 = ', collisions_1

print 'collisions_2 = ', collisions_2

collisions_1 = [i for i in chain.from_iterable(collisions_1)]
print  'collisions in calcite 1 = ', collisions_1 

collisions_2 = [i for i in chain.from_iterable(collisions_2)]
print 'collisions in calcite 2 = ', collisions_2 

output_array_1 = np.vstack((collisions_1))
np.savetxt('collisions_1.dat', output_array_1, header="P(GPa) \t   H per F unit (a.u)", fmt="%0.13f")

output_array_2 = np.vstack((collisions_2))
np.savetxt('collisions_2.dat', output_array_2, header="P(GPa) \t   H per F unit (a.u)", fmt="%0.13f")

Intersection = max(collisions_1, key=lambda item:item[1])
print ' Intersection = ', Intersection


# Plotting and fitting to the exact expression of Delta_H:
#********* Exact expression of Delta_H:

fig = plt.figure()

EnergyCI, VolumeCI, PressureCI, EnthalpyCI  = np.loadtxt('./E_V_P_H__C_I.dat', skiprows = 1).T

print 'PressureCI[0] = ', PressureCI[0]
print 'PressureCI[-1] = ', PressureCI[-1]

Energy14, Volume14, Pressure14, Enthalpy14  = np.loadtxt('./E_V_P_H__14.dat', skiprows = 1).T

print 'Pressure14[0] = ', Pressure14[0]
print 'Pressure14[-1] = ', Pressure14[-1]

#xp_C_I = np.linspace(PressureCI[-1], PressureCI[0], 100)
xp_C_I = np.linspace(PressureCI[0], PressureCI[-1], 100)
#xp_14 = np.linspace(Pressure14[-1], Pressure14[0], 100)
xp_14 = np.linspace(Pressure14[0], Pressure14[-1], 100)

# Plotting the fitting curves:

# This does not work (yields a wrong fit, because you plot plot H(V) versus some completely uncorrelated pressures xp_C_I:
#p2, = plt.plot(xp_C_I, H(V_C_I_lin, *popt_C_I), color='grey', label='H fit Calcite I' )

# This works but it does not play with all the linspace in plotting 
#p2, = plt.plot(PressureCI, H(V_C_I, *popt_C_I), color='grey', label='H fit Calcite I' )

# IOBE suggestion:
p2, = plt.plot(P(V_C_I_lin, *popt_C_I), H(V_C_I_lin, *popt_C_I), color='black', label='H fit Data' )

# This does not work (yields a wrong fit, because you plot plot H(V) versus some completely uncorrelated pressures xp_C_I:
#p6, = plt.plot(xp_14, H(V_14_lin, *popt_14), 'b', label='H fit Calcite II')

# This works but it does not play with all the linspace in plotting
#p6, = plt.plot(Pressure14, H(V_14, *popt_14), 'b', label='H fit Calcite II')

# IOBE suggestion:
p6, = plt.plot(P(V_14_lin, *popt_14), H(V_14_lin, *popt_14), color='blue', label='H fit Data' )

print 'V_C_I_lin = ', V_C_I_lin
print ' xp_C_I = ', xp_C_I
print ' H(V_C_I_lin, *popt_C_I) = ', H(V_C_I_lin, *popt_C_I)

# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

fontP = FontProperties()
fontP.set_size('small')

plt.legend((p1, p2, p5, p6), ("Calcite I", 'H fit Calcite I', 'Calcite II', 'H fit Calcite II'), prop=fontP)

plt.xlabel('P / Formula unit (GPa)')
plt.ylabel('H / Formula unit (a.u.)')
plt.suptitle("B3LYP, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
plt.ticklabel_format(useOffset=False)
ax = fig.add_subplot(111)
ax.annotate('Intersection\nP= %g GPa\nH = %g a.u.' %(Intersection[0], Intersection[1]), xy=(Intersection[0], Intersection[1]), xytext=(Intersection[0]+2.7767, Intersection[1]-162.27),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='purple'),
            )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_exact_expression_with_intersection.pdf', bbox_inches='tight')








plt.show()

#sys.exit()




