import matplotlib
matplotlib.use('TKagg')
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
from matplotlib.path import Path
from dirac_sheet import dirac_sheet


def main():
	plt.ion()
	
	fig = plt.figure(num=None, figsize=(11, 6), dpi=96, facecolor='w', edgecolor='k')
	
	#imput parameters for Dirac_sheet
	dx = .25
	dt = 0.1
	Ngrid = 1802
	X_offset=112.5
	Nsteps=2200
	z_lim=2
	filename='Dirac_run_out'

	plt.clf()

	tic = time.time()
	
	myDirac=dirac_sheet(0,(901,601),dt,dx,X_offset,0)
	myDirac.set_p(.2,np.pi*0/180)
	myDirac.massorV=False
	x, y = myDirac.get_pos_mat()
	

	NoPropMat = np.zeros(x.shape,dtype=np.uint8)
	AbsMat = np.zeros(x.shape)
	DriveMat=np.zeros(x.shape)
	
	#define part of injector shape
	poly_verts = np.array([[0,37],[70,12.5],[70,1001],[0,1001],[0,37]])
	poly_verts2 = np.array([[0,37],[70,12.5],[70,13.5],[0,38],[0,37]])
	
	NoPropMat[inPolygon(poly_verts2,x,y)]=1
	AbsMat[inPolygon(poly_verts,x,y)]=1
	#mirror across x-axis and make CCW
	poly_verts[:,1]*=-1
	poly_verts[:,1]-=.125
	poly_verts2[:,1]*=-1
	poly_verts2[:,1]-=.125
	
	poly_verts[:,:]=poly_verts[::-1,:]
	poly_verts2[:,:]=poly_verts2[::-1,:]
	
	NoPropMat[inPolygon(poly_verts2,x,y)]=1
	AbsMat[inPolygon(poly_verts,x,y)]=1
	NoPropMat[((x<140.26)&(x>140)&(y<12.5)&(y>-12.5))]=1
	
	AbsMat[x>140]=.99
	AbsMat[(x>70)&(x<141)&(y>40)]=.99
	AbsMat[(x>70)&(x<141)&(y<-40)]=.99
	
	DriveMat[(x>0)&(x<1)&(y>-36)&(y<36)]=1
	
	myDirac.set_No_prop_mat(NoPropMat)
	myDirac.set_Absorb_mat(AbsMat)
	myDirac.set_Drive_mat(DriveMat)
	
	
	for i in xrange(0,Nsteps):
		#timestep
		myDirac.time_step()
		toc = time.time()
		print i, toc-tic
		tic = time.time()
		
		
		if np.mod(i,10)==0:
			#print colormap every 10 step
			plt.clf()
			
			
			#z limits
			z_min, z_max = -z_lim, z_lim #np.min(np.real(myDirac.u1)), np.max(np.real(myDirac.u1))


			
			plt.subplot(1,2, 1)
			plt.imshow(np.real(myDirac.u1)-0*np.real(myDirac.u10*np.exp(-1j*myDirac.p0*myDirac.t)), cmap='RdBu', vmin=z_min, vmax=z_max,
					   extent=[x.min(), x.max(), y.min(), y.max()],
					   interpolation='nearest', origin='lower')
			plt.title('u1')
			ax = plt.gca()
			ax.set_xlim([0,200])
			#ax.set_ylim([-25,25])
			plt.colorbar()
			plt.subplot(1,2,2)
			plt.imshow(np.real(myDirac.u10*np.exp(-1j*myDirac.p0*myDirac.t)), cmap='RdBu', vmin=z_min, vmax=z_max,
					   extent=[x.min(), x.max(), y.min(), y.max()],
					   interpolation='nearest', origin='lower')
			plt.title('v1')
			ax = plt.gca()
			ax.set_xlim([0,200])
			#ax.set_ylim([-25,25])
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