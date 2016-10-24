import matplotlib
matplotlib.use('TKagg')
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
from matplotlib.path import Path
import dirac_script as ds
import multiprocessing
import os


def main():
	directory = "E:/pyDirac/collimator_output_full_absorb/"
	#
	
	plt.ion()
	fig = plt.figure(figsize=(12, 6), dpi=96, facecolor='w', edgecolor='k')

	plt.clf()

	ax = plt.plot()
	
	for wave_vec in np.linspace(.15,1.5,10,endpoint=True):
		for theta in np.linspace(-30.0,30.0,21,endpoint=True):
			#theta=0;
			tic = time.time()
			
			data=np.load(directory+'collimator_%03d_deg_%03dEm3_wave_vec.npz' %(theta,wave_vec*100))
	
			u1=data['u1']
			#u2=data['u2']
			v1=data['v1']
			#v2=data['v2']
			plt.clf()

			#z limits
			z_min, z_max = -3, 3 #np.min(np.real(myDirac.u1)), np.max(np.real(myDirac.u1))


			#plt.ioff()
			plt.subplot(1,2, 1)
			plt.imshow(np.real(u1).T, cmap='RdBu', vmin=z_min, vmax=z_max,
					   #extent=[x.min(), x.max(), y.min(), y.max()],
					   interpolation='nearest', origin='lower')
			plt.title('u1')
			plt.colorbar()
			plt.subplot(1,2,2)
			plt.imshow(np.real(v1).T, cmap='RdBu', vmin=z_min, vmax=z_max,
					   #extent=[x.min(), x.max(), y.min(), y.max()],
					   interpolation='nearest', origin='lower')
			plt.title('v1')
			plt.colorbar()
			
			fig.canvas.draw()
			#plt.draw()
		
			
			
				
			toc = time.time()
				
			print toc-tic
	
			
			
			
if __name__ == '__main__':
	main()