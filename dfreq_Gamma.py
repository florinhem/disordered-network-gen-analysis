
"""
w=ck in vaccum
w=|k|c/n_eff in homogenous medium
w_G=|k_G|c/n_eff for Gamma point
k_G=b(h,k,l)
b=2 pi/a
a=l_s*d for supercell_length times rod length
n_eff=n_1*phi_1+n_2*phi_2

=> wd/2 pi c = 1/(l_s*n_eff) *(h^2+k^2+l^2)
"""
import numpy as np

network_type="ctn"
n_air=1
n_cylinder=3.6
fill_fraction=n_air/(n_air+n_cylinder)
n_eff=n_cylinder*fill_fraction+n_air*(1-fill_fraction)

print(network_type)
print(fill_fraction)
print(n_eff)

if network_type=="dia":
    supercell_length=4/np.sqrt(3)
elif network_type=="srs":
    supercell_length=4/np.sqrt(2)
elif network_type=="srd":
    supercell_length=2.3075
elif network_type=="pto":
    supercell_length=2.8284
elif network_type=="ctn":
    supercell_length=3.7033
else:
    print("unknown network type")

n=3
r=range(-n,n)
k_Gamma_array=[(h,k,l) for h in r for k in r for l in r]

print()
print(k_Gamma_array)
print()

norm_freq=[1/(supercell_length*n_eff)*np.linalg.norm(k_Gamma) for k_Gamma in k_Gamma_array]

print()
print(norm_freq)
print()

sort_norm_freq=np.sort(np.round(norm_freq,decimals=3))

print()
print(sort_norm_freq)
print()

norm_freq_dict=[
    (np.round(sort_norm_freq[i+1]-sort_norm_freq[i],decimals=3),
     np.sum(sort_norm_freq==sort_norm_freq[i])) 
     for i in range(len(sort_norm_freq)-1) 
     if sort_norm_freq[i+1]!=sort_norm_freq[i]]

print()
print(norm_freq_dict)
print()

spacing_norm_freq_dict=[
    ( np.round(sort_norm_freq[i+1]-sort_norm_freq[i],decimals=3),
     (np.sum(sort_norm_freq==sort_norm_freq[i+1])+np.sum(sort_norm_freq==sort_norm_freq[i]))/2) 
     for i in range(len(sort_norm_freq)-1) 
     if sort_norm_freq[i+1]!=sort_norm_freq[i]]

print()
print(spacing_norm_freq_dict)
print()

distance=[delta_omega/nr_bands for (delta_omega, nr_bands) in spacing_norm_freq_dict]

print()
print(distance)
print()

#distance_norm_freq=[sort_norm_freq[i+1]-sort_norm_freq[i] for i in range(len(sort_norm_freq)-1)]
#print(distance_norm_freq)