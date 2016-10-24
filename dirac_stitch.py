import numpy as np
import time
from dirac_sheet import dirac_sheet
pi=np.pi






class dirac_stitch:

	def __init__(self,sheet1,sheet2)
		self.sheet1=sheet1
		self.sheet2=sheet2
		x1, y1 = sheet1.get_pos_mat()
		x2, y2 = sheet2.get_pos_mat()
		
		self.D_x=sheet1.D_x
		
		self.nXoffset=np.round((x2[0,0]-x1[0,0])/self.D_x)
		self.nYoffset=np.round((y2[0,0]-y1[0,0])/self.D_x)
		
		