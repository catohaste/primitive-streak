from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

number_of_cells = 750

fig = plt.figure()
ax =  plt.axes(xlim=(0,number_of_cells),ylim=(0,7),xlabel='Cell Location')

line1_cells = np.arange(0,number_of_cells/2,1)
line2_cells = np.arange(number_of_cells/2,number_of_cells,1)
line1_i = np.zeros((int(number_of_cells/2),1),dtype=float)
line2_i = np.zeros((int(number_of_cells/2),1),dtype=float)

line_cells = np.arange(0,number_of_cells,1)
line_i = np.zeros((number_of_cells,1),dtype=float)

a = 0.000007

counter = 0
for cell in line1_cells:
    line1_i[counter] = (6/number_of_cells * cell) + 3
    counter = counter+1

counter = 0
for cell in line2_cells:
    line2_i[counter] = (-6/number_of_cells * cell) + 9
    counter = counter+1
    
counter = 0
for cell in line_cells:
    line_i[counter] = -a*np.power(cell-(number_of_cells-1)/2,2) + 2
    counter = counter+1

ax.plot(line1_cells,line1_i,'r-',markersize = 0.5)
ax.plot(line2_cells,line2_i,'b-',markersize = 0.5)
ax.plot(line_cells,line_i,'k-',markersize = 0.5)
ax

plt.show()