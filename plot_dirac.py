import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
from matplotlib.path import Path
from dirac_sheet import dirac_sheet


def main():
	fig = plt.figure()

	# make these smaller to increase the resolution

	dx = .25
	dt = 0.1
	Ngrid = int(1802)
	X_offset=112.5
	# generate 2 2d grids for the x & y bounds
	#y, x = np.ogrid[slice(-5, 5 + dy, dy),	
	#                slice(-3, 3 + dx, dx)]


	#y, x = np.mgrid[slice((0-round(Ngrid/2))*dy,(Ngrid-round(Ngrid/2))*dy,dy),
	#				slice((0-round(Ngrid/2))*dx+112.25,(Ngrid-round(Ngrid/2))*dx+112.25,dx)]

	plt.ion()

	plt.clf()


	ax = plt.plot()

	
	# tic = time.time()
	
	myDirac=dirac_sheet(0,901,dt,dx,X_offset,0)
	myDirac.set_p(.5,np.pi/4.0)
	
	# x, y = myDirac.get_pos_mat()
	

	# NoPropMat = np.zeros(x.shape,dtype=np.uint8)
	
	# print myDirac.p0
	
	#define part of injector shape
	#poly_verts = np.array([[0,37],[70,12.5],[70,1001],[0,1001],[0,37]])
	
	#NoPropMat[inPolygon(poly_verts,x,y)]=1
	# mirror across x-axis and make CCW
	# poly_verts[:,1]*=-1
	# poly_verts[:,1]-=.125
	
	# poly_verts[:,:]=poly_verts[::-1,:]
	
	# NoPropMat[inPolygon(poly_verts,x,y)]=1
	# NoPropMat[((x<140.26)&(x>140)&(y>12.5))]=1
	# NoPropMat[((x<140.26)&(x>140)&(y<-12.5))]=1
	
	# AbsMat = np.zeros(x.shape)
	# AbsMat[x>205]=.99
	# AbsMat[(x>70)&(x<141)&(y>40)]=.99
	# AbsMat[(x>70)&(x<141)&(y<-40)]=.99
	
	# DriveMat=np.zeros(x.shape)
	# DriveMat[(x>0)&(x<1)&(y>-36)&(y<36)]=1
	
	# myDirac.set_No_prop_mat(NoPropMat)
	# myDirac.set_Absorb_mat(AbsMat)
	# myDirac.set_Drive_mat(DriveMat)
	
	# data=np.load('file.npz')
	# u1=data['u1']
	# u2=data['u2']
	# v1=data['v1']
	# v2=data['v2']
	# plt.clf()

	#z limits
	z_min, z_max = -0.3, 0.3 #np.min(np.real(myDirac.u1)), np.max(np.real(myDirac.u1))

	print np.max(np.real((myDirac.u10*np.conj(myDirac.v10))+(myDirac.v10*np.conj(myDirac.u10))).T)
	print np.max(np.imag((myDirac.u10*np.conj(myDirac.v10))-(myDirac.v10*np.conj(myDirac.u10))).T)
	print myDirac.theta
	
	plt.imshow(np.real(myDirac.u10),cmap='RdBu', vmin=z_min, vmax=z_max,
	#plt.imshow(np.real(v1*np.conj(v1)).T+np.real(u2*np.conj(u2)).T+0*np.real(v1*np.conj(v1)).T+0*np.real(v2*np.conj(v2)).T, cmap='hot', vmin=z_min, vmax=z_max,
	#plt.imshow(np.imag((myDirac.u10*np.conj(myDirac.v10))-(myDirac.v10*np.conj(myDirac.u10))).T, cmap='RdBu', vmin=z_min, vmax=z_max,
	#plt.imshow(np.real(((u1+u2)*np.conj(v1+v2))+((v1+v2)*np.conj(u1+u2))).T, cmap='RdBu', vmin=z_min, vmax=z_max,
	#plt.imshow(np.imag(np.gradient(u2)[1]*(np.conj(u2))-u2*np.gradient(np.conj(u2))[1]).T, cmap='RdBu', vmin=z_min, vmax=z_max,
			   #extent=[x.min(), x.max(), y.min(), y.max()],
			   interpolation='nearest', origin='lower')
	#plt.quiver(np.imag(np.gradient(u2)[0]*(np.conj(u2))-u2*np.gradient(np.conj(u2))[0]).T,np.imag(np.gradient(u2)[1]*(np.conj(u2))-u2*np.gradient(np.conj(u2))[1]).T)
	#plt.title('image (interp. nearest)')	
	#plt.colorbar()

	
	fig.canvas.draw()
	
	
	plt.ioff()
	# This line was moved up <----
	plt.draw()
	plt.show()
		
		
def inPolygon(poly_verts,x,y):

	ss=x.shape
	
	
	x1, y1 = x.flatten(), y.flatten()

	points = np.vstack((x1,y1)).T

	path = Path(poly_verts)
	
	grid = path.contains_points(points, radius=.01)
	grid = grid.reshape(ss)
	
	return grid
		
if __name__ == '__main__':
	main()