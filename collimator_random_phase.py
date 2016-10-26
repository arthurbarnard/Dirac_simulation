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
	directory = "E:/pyDirac/collimator_output_full_absorb_2/"
	
	# plt.ion()
	# fig = plt.figure(figsize=(12, 6), dpi=96, facecolor='w', edgecolor='k')

	# plt.clf()

	ax = plt.plot()
	
	theta=-30
	wave_vec=.15
	data=np.load(directory+'collimator_%03d_deg_%03dEm3_wave_vec.npz' %(theta,np.round(wave_vec*100)))
	u1=data['u1']
	u1s=np.zeros((u1.shape[0],u1.shape[1],21),dtype=np.complex)
	v1s=np.zeros((u1.shape[0],u1.shape[1],21),dtype=np.complex)
	
	for wave_vec in np.linspace(.15,1.95,13,endpoint=True):
		#wave_vec=1.05
		for n_index, theta in enumerate(np.linspace(-30.0,30.0,21,endpoint=True)):
			#theta=0;
			tic = time.time()
			data=np.load(directory+'collimator_%03d_deg_%03dEm3_wave_vec.npz' %(theta,np.round(wave_vec*100)))

			u1=data['u1']
			v1=data['v1']
			u1s[:,:,n_index]=u1
			v1s[:,:,n_index]=v1

		dtheta=15
			
		psi=np.zeros(u1.shape)
		jx=np.zeros(u1.shape)
		jx2=np.zeros(u1.shape)
		
		fname=directory+'collimator2_%03dEm3_wave_vec' %(np.round(wave_vec*100))
		print fname
		if not os.path.isfile(fname+'.npz'):
			for k in xrange(1000):
				u1[:]=0
				v1[:]=0
				pos=np.random.rand()*72-36.0
				chirality=np.round(np.random.rand())==1
				for n_index, theta in enumerate(np.linspace(-30.0,30.0,21,endpoint=True)):
					#u1+=u1s[:,:,n_index]*np.exp(-(theta-theta0)**2/2/dtheta**2)
					#v1+=v1s[:,:,n_index]*np.exp(-(theta-theta0)**2/2/dtheta**2)
					amp=np.random.rand()
				
					u1+=u1s[:,:,n_index]*amp*np.exp(1j*pos*np.sin(theta/180*np.pi))
					if chirality:
						v1+=v1s[:,:,n_index]*amp*np.exp(1j*pos*np.sin(theta/180*np.pi))
					else:
						v1+=v1s[::-1,:,20-n_index]*amp*np.exp(1j*pos*np.sin(theta/180*np.pi))
					
				psi+=np.real(u1*np.conj(u1))
				u1x=np.gradient(u1)[1]
				jx-=np.imag(u1*np.conj(u1x)-np.conj(u1)*u1x)
				jx2+=np.real(u1*np.conj(v1)+np.conj(u1)*v1)
				# plt.clf()

				# #z limits
				# z_min, z_max = -15, 15 #np.min(np.real(myDirac.u1)), np.max(np.real(myDirac.u1))


				# plt.ioff()
				# plt.subplot(1,2, 1)
				# plt.imshow(np.real(u1s[:,:,10]).T*4, cmap='RdBu', vmin=z_min, vmax=z_max,
						   # #extent=[x.min(), x.max(), y.min(), y.max()],
						   # interpolation='nearest', origin='lower')
				# #plt.title('u1')
				# plt.axis('off')
				# #plt.colorbar()
				# ax=plt.gca()
				# ax.set_ylim([0, 800])
				# plt.subplot(1,2,2)
				# plt.imshow(np.real(jx2).T, cmap='RdBu', #vmin=z_min, vmax=z_max,
						   # #extent=[x.min(), x.max(), y.min(), y.max()],
						   # interpolation='nearest', origin='lower')
				# #plt.title('v1')
				# plt.axis('off')
				# #plt.colorbar()
				# ax=plt.gca()
				# ax.set_ylim([0, 800])
				
				# fig.canvas.draw()
				#plt.draw()
			
				
				
					
				toc = time.time()
					
				print toc-tic
			
			
		
				
			np.savez_compressed(fname+'.npz', psi=psi,jx=jx,jx2=jx2)	

		
if __name__ == '__main__':
	main()