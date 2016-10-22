import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
from matplotlib.path import Path
from dirac_sheet import dirac_sheet


def main():
	fig = plt.figure(num=None, figsize=(11, 6), dpi=96, facecolor='w', edgecolor='k')
	#small comment added
	#imput parameters for Dirac_sheet
	dx = .25
	dt = 0.1
	Ngrid = 1802
	X_offset=112.5
	Nsteps=2200
	z_lim=2
	filename='Dirac_run_out'

	plt.ion()
	plt.clf()

	tic = time.time()
	
	myDirac=dirac_sheet(0,901,dt,dx,X_offset,0)
	myDirac.set_p(.1,np.pi*0)
	
	x, y = myDirac.get_pos_mat()
	

	NoPropMat = np.zeros(x.shape,dtype=np.uint8)
	
	#define part of injector shape
	poly_verts = np.array([[0,37],[70,12.5],[70,1001],[0,1001],[0,37]])
	
	NoPropMat[inPolygon(poly_verts,x,y)]=1
	
	# mirror across x-axis and make CCW
	poly_verts[:,1]*=-1
	poly_verts[:,1]-=.125
	poly_verts[:,:]=poly_verts[::-1,:]
	NoPropMat[inPolygon(poly_verts,x,y)]=1
	
	#define second aperture
	NoPropMat[((x<140.26)&(x>140)&(y>12.5))]=1
	NoPropMat[((x<140.26)&(x>140)&(y<-12.5))]=1
	
	#Define absorption matrix
	AbsMat = np.zeros(x.shape)
	AbsMat[x>205]=.99
	AbsMat[(x>70)&(x<141)&(y>40)]=.99
	AbsMat[(x>70)&(x<141)&(y<-40)]=.99
	
	DriveMat=np.zeros(x.shape)
	DriveMat[(x>0)&(x<1)&(y>-36)&(y<36)]=1
	
	myDirac.set_No_prop_mat(NoPropMat)
	myDirac.set_Absorb_mat(AbsMat)
	myDirac.set_Drive_mat(DriveMat)
	
	
	for i in xrange(0,Nsteps):
		#timestep
		myDirac.time_step()
		toc = time.time()
		print i, toc-tic
		
		
		if np.mod(i,10)==0:
			#print colormap every 10 step
			plt.clf()
			
			
			#z limits
			z_min, z_max = -z_lim, z_lim #np.min(np.real(myDirac.u1)), np.max(np.real(myDirac.u1))


			
			plt.subplot(1,2, 1)
			plt.imshow(np.real(myDirac.u1), cmap='RdBu', vmin=z_min, vmax=z_max,
					   extent=[x.min(), x.max(), y.min(), y.max()],
					   interpolation='nearest', origin='lower')
			plt.title('u1')
			ax = plt.gca()
			ax.set_xlim([0,200])
			plt.colorbar()
			plt.subplot(1,2,2)
			plt.imshow(np.real(myDirac.v1), cmap='RdBu', vmin=z_min, vmax=z_max,
					   extent=[x.min(), x.max(), y.min(), y.max()],
					   interpolation='nearest', origin='lower')
			plt.title('v1')
			ax = plt.gca()
			ax.set_xlim([0,200])
			plt.colorbar()
			
			fig.canvas.draw()
			
	
	#save key variables from simulation for later analysis
	np.savez_compressed(filename+'.npz', u1=myDirac.u1,u2=myDirac.u2,v1=myDirac.v1,v2=myDirac.v2,Ngrid=myDirac.Ngrid,t=myDirac.t,D_x=myDirac.D_x)
	plt.ioff()
	plt.draw()
	plt.show()
		
		
def inPolygon(poly_verts,x,y):
	#outputs logical array for which x,y are inside poly_verts
	ss=x.shape
	
	x1, y1 = x.flatten(), y.flatten()

	points = np.vstack((x1,y1)).T

	path = Path(poly_verts)
	
	grid = path.contains_points(points, radius=.01)
	grid = grid.reshape(ss)
	
	return grid
		
if __name__ == '__main__':
	main()