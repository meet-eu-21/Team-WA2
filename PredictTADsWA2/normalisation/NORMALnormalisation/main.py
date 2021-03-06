#python 3
#2020
#CC-By-SA
#Carron Leopold, Julien Mozziconacci

import sys
import numpy as np
from scipy.spatial import ConvexHull
from scipy import sparse
import pandas as pd
from sklearn.manifold import MDS #If you want the scikit learn mds
import HiCtoolbox
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-dir', '--inputdir', help='directory with input HiC contact matrices for chromosomes.')
parser.add_argument('-chr', '--chromosomenr', help='number of the chromosome to compute normalised matrix.')
#parser.add_argument('-path', '--path', help='path to the RAWobserved files.')

#args = parser.parse_args()
args = vars(parser.parse_args())
chromosome_nr = "chr%s" % args['chromosomenr']
path = args['inputdir']
#print("Normalising chromosome %s matrix" %chromosome_nr)
#path = "/home/eu_plus_3/HiC/"

R=100000
NbmaxEpi=15 #Epi states go from 0 to 15
alpha=0.227
selectedmark=1 #index of the selected mark
res = path.split("/")[-2].split("_")[0]
HiCfilename=path+chromosome_nr+'_'+res+'.RAWobserved'
EpiGfilename='E116_15_coreMarks_dense'

#Build matrix
A=np.loadtxt(HiCfilename)
A=np.int_(A)
#print('Input data shape : ',np.shape(A))
A=np.concatenate((A,np.transpose(np.array([A[:,1],A[:,0],A[:,2]]))), axis=0)#build array at pb resolution
A = sparse.coo_matrix( (A[:,2], (A[:,0],A[:,1])))
binned_map=HiCtoolbox.bin2d(A,R,R) #!become csr sparse array
#print(binned_map)
LENTEST=np.shape(A)[0]
#print('Input at the good resolution : ',np.shape(binned_map))

del A #keep space

#Build color annotation at desired resolution
color=pd.read_csv(EpiGfilename,delimiter='\t',header=None,names=[1,2,3,4])
color=color[color[1]==chromosome_nr]#take only chr of interest
#print(color)
number=color[4].max() #number of color in the file
color_vec=np.zeros((LENTEST,number+1)) #build array at pb resolution LENchr * number of color
i=0
while i<np.shape(color)[0]:
	color_vec[color[2].iloc[i]:color[3].iloc[i],color[4].iloc[i]]=1
	i+=1

color_bins=HiCtoolbox.bin2d(color_vec,R,1)
color_bins=color_bins/np.amax(color_bins)

#print('Bp cover by this mark, has to be >0 :',np.sum(color_bins[:,selectedmark]) )

# FILTER
#print("before filtering : ",np.shape(binned_map))
sumHicmat=np.sum(binned_map,0) 
mini = np.mean(sumHicmat)-np.std(sumHicmat)*1.5 #min value of filtering
maxi = np.mean(sumHicmat)+np.std(sumHicmat)*1.5 #max value of filtering
binsaved=np.where(np.logical_and(mini < sumHicmat,sumHicmat < maxi)) #coord of bin to save
filtered_map=binned_map[binsaved[1],:] #save on row
filtered_map=filtered_map[:,binsaved[1]] #save on col
#print("after filtering :??",np.shape(filtered_map))#,np.shape(color_vecseg))

color2=color_bins[binsaved[1]] #filter the epi by removed bin in HiC
color2=color2[:,selectedmark] #now color2 is 1D
color2=np.float64(color2.todense()) #type issue

#3D
#print('3D')#Here : sparse int64
contact_map=HiCtoolbox.SCN(filtered_map.copy()) 
contact_map=np.asarray(contact_map)**alpha #now we are not sparse at all

df = pd.DataFrame(contact_map)
# dodawanie trzech kolumn potrzebnych do orginalnego TopDom, ale mojej wersji juz nie
# df.insert(0, chromosome_nr, chromosome_nr)
# df.insert(1, 'pos', np.arange(0,len(df.index)*100000,100000))
# df.insert(2, 'end', np.arange(100000,(len(df.index)+1)*100000,100000))
df.to_csv('normalised_matrices/contact_map_%s.tsv' %chromosome_nr, sep='\t', header=None, index=False)
print('normalised_matrices/contact_map_%s.tsv' %chromosome_nr)

#print(contact_map.shape)
plt.imshow(contact_map)
plt.show()
# dist_matrix = HiCtoolbox.fastFloyd(1/contact_map) #shortest path on the matrix
# dist_matrix=dist_matrix-np.diag(np.diag(dist_matrix))#remove the diagonal
# dist_matrix=(dist_matrix+np.transpose(dist_matrix))/2; #just to be sure that the matrix is symetric, not really usefull in theory


# #MDS
# embedding = MDS(n_components=3)#LOAD the MDS #With scikit-learn mds
# XYZ = embedding.fit_transform(dist_matrix) #Make the transform
# XYZ=np.float64(XYZ)
# XYZ,E=HiCtoolbox.sammon(dist_matrix, 3)#with the one from tom j pollard




# print("Output shape : ",np.shape(XYZ),np.shape(color2))

# #point rescale
# hull=ConvexHull(XYZ)
# scale=100/hull.area**(1/3)
# XYZ=XYZ*scale
# print(XYZ)
# HiCtoolbox.writePDB('3Dcolors_'+str(alpha)+'.pdb',XYZ,color2)
