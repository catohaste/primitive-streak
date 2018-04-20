import matplotlib.pyplot as plt
import numpy as np

number_of_cells = 751

# fixed model parameters
inhibitor_m = 0.000022
inhibitor_c = 1

diffusion_inhibitor = 0

# fixed values of the 'variables'
inhibitor_shift = 5 
propagator = 467

# initiate parameters to vary
#SHORT decay_inhibitor = [0.002] # per second
#LONG 
decay_inhibitor = np.concatenate( ( np.arange(0.001,0.01,0.0005), np.arange(0.01,0.02,0.001), np.arange(0.02,0.1,0.005),np.arange(0.1,0.4,0.05) ) )
#SHORT saturation_i = [0.01] # micro meters per second
#LONG 
saturation_i = np.concatenate( (  np.arange(0.01,0.1,0.005), np.arange(0.1,1,0.05), np.arange(1,2,0.1) ) )
#SHORT k_i = np.arange(100,1000,50)  # nano molar
#LONG 
k_i = np.concatenate((np.arange(100,1000,50),np.arange(1000,2000,100),np.arange(2000,10000,1000)))
#SHORT a_ii = np.arange(10,100,5) # no units
#LONG
a_ii = np.concatenate((np.arange(10,100,5),np.arange(100,200,10),np.arange(200,1000,50))) 


number_of_points = len(decay_inhibitor) * len(saturation_i) * len(k_i) * len(a_ii)

print(number_of_points)

dV = np.zeros((len(decay_inhibitor),len(saturation_i),len(k_i),len(a_ii)),dtype=np.float32)

counter = 0
for idx_decay,decay in enumerate(decay_inhibitor):
    for idx_sat,sat in enumerate(saturation_i):
        for idx_k,k in enumerate(k_i):
            for idx_a,a in enumerate(a_ii):
                dV[idx_decay,idx_sat,idx_k,idx_a] = np.absolute( ( (sat*( np.power(a * inhibitor_shift,4) + np.power(propagator,4) ) ) / ( k + np.power(a*inhibitor_shift,4) + np.power(propagator,4) ) ) - (decay * inhibitor_shift) )
                # if np.absolute(dV[idx_decay,idx_sat,idx_k,idx_a]) < 0.0000000000001:
#                     counter = counter + 1
#                     print((idx_decay,idx_sat,idx_k,idx_a))
#
# print(counter)


ind = np.unravel_index(np.argmin(np.absolute(dV), axis=None), dV.shape)
print( ind , dV[ind] )

print('decay =',decay_inhibitor[ind[0]])
print('saturation =',saturation_i[ind[1]])
print('K =',k_i[ind[2]])
print('a_ii =',a_ii[ind[3]])



#print(dV[0,0,0,0])
#print(decay_inhibitor[2],saturation_i[0])
#print('ratio =',saturation_i[0] / decay_inhibitor[2])

# print(dV[48,41,8,10])
# print(decay_inhibitor[48],saturation_i[41])
# print('ratio =',saturation_i[41] / decay_inhibitor[48] )

plot_data = dV[:,0,:,:]


fig, ax = plt.subplots()
im = plt.imshow(plot_data,cmap='hot',interpolation='nearest')
plt.colorbar()

ax.set_xticks(range(len(a_ii)))
ax.set_yticks(range(len(k_i)))

ax.set_xticklabels(a_ii,rotation='vertical')
ax.set_yticklabels(k_i)


ax.set_xlabel('')
ax.set_ylabel('')

fig.tight_layout()
plt.show()


# print(dV[2,0,8,10])
# print(decay_inhibitor[2],saturation_i[0])
# print('ratio =',saturation_i[0] / decay_inhibitor[2])
#
# print(dV[48,41,8,10])
# print(decay_inhibitor[48],saturation_i[41])
# print('ratio =',saturation_i[41] / decay_inhibitor[48] )
#
# print(len(decay_inhibitor))
# print(len(saturation_i))
#
# plot_data_decay_sat = dV[:,:,8,10]
#
#
# fig_decay_sat, ax_decay_sat = plt.subplots()
# im_decay_sat = plt.imshow(plot_data_decay_sat,cmap='hot',interpolation='nearest')
# plt.colorbar()
#
# ax_decay_sat.set_xticks(range(len(saturation_i)))
# ax_decay_sat.set_yticks(range(len(decay_inhibitor)))
#
# ax_decay_sat.set_xticklabels(np.round(saturation_i,decimals = 2),rotation='vertical')
# ax_decay_sat.set_yticklabels(np.round(decay_inhibitor,decimals = 3))
#
#
# ax_decay_sat.set_xlabel('saturation')
# ax_decay_sat.set_ylabel('decay')
#
# fig_decay_sat.tight_layout()
# plt.show()




