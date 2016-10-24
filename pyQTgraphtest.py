import pyqtgraph as pg
import matplotlib.pyplot as plt
import numpy as np

def main():
	y,x = np.ogrid[slice(-10,10,100),slice(-10,10,100)]
	z=np.exp(1j*x/2/np.pi)
	pg.image(z)
	
	
if __name__ == '__main__':
	main()