import matplotlib.pyplot as plt
import numpy as np

number_of_cells = 751

fig = plt.figure()
ax =  plt.axes(xlim=(0,number_of_cells),ylim=(0,7),xlabel='Cell Location')

line_cells = np.arange(0,number_of_cells,1)
line_i = np.zeros((number_of_cells,1),dtype=float)

m = 0.000022
c = 11
    
counter = 0
for cell in line_cells:
    line_i[counter] = c - m*np.power(cell-(number_of_cells-1)/2,2)
    counter = counter+1

ax.plot(line_cells,line_i - 5,'k-')
ax

plt.show()