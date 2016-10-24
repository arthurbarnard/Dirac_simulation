import dirac_script_absorb_everywhere as ds
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
from matplotlib.path import Path
import multiprocessing
import os


def main():

	directory = "E:/pyDirac/collimator_output_full_absorb_2/"
	if not os.path.exists(directory):
		os.makedirs(directory)
	
	jobs=[]
	n=0
	for wave_vec in np.linspace(0.15,1.5,10,endpoint=True):		
		for theta in np.linspace(-30.0,30.0,21,endpoint=True):

			p = multiprocessing.Process(target=ds.run_Dirac, args=(wave_vec,theta*np.pi/180.0,directory+'collimator_%03d_deg_%03dEm3_wave_vec' %(theta,wave_vec*100)))
			jobs.append(p)
			p.start()
			n+=1
			print 'collimator_%03d_deg_%03dEm3_wave_vec' %(theta,wave_vec*100)
			if np.mod(n,7)==0:
				p.join()
				
if __name__ == '__main__':
	main()