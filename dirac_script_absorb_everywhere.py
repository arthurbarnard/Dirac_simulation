from dirac_run import inPolygon
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
from matplotlib.path import Path
from dirac_sheet import dirac_sheet


def main():
	
	run_Dirac(.1,0,'file_out')

	
def run_Dirac(p,theta,fname):

	dx = .25
	dt = 0.1
	Ngrid = int(1802)
	X_offset=112.5

	tic = time.time()
	
	myDirac=dirac_sheet(0,(901,601),dt,dx,X_offset,0)
	myDirac.set_p(p,theta)
	
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
	

	
	for i in xrange(0,2200):
		
		myDirac.time_step()
		toc = time.time()
		print toc-tic
		
		tic = time.time()
		
			
	np.savez_compressed(fname+'.npz', u1=myDirac.u1,u2=myDirac.u2,v1=myDirac.v1,v2=myDirac.v2,Ngrid=(myDirac.Ngrid_x,myDirac.Ngrid_y),t=myDirac.t,D_x=myDirac.D_x)	
		
		
if __name__ == '__main__':
	main()