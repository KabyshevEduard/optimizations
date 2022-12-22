import numpy as np

def f(x1, x2):
	return x1**2 + x2**2 - 1.2*x1*x2

def rosen(x1, x2):
	return (1-x1)**2 + 100 * (x2 - x1**2)**2
 
def dnf_dxn(X, n, f):
	h = 0.01
	df_dx1 = (f(X[0]+h/2, X[1]) - f(X[0]-h/2, X[1]))/h
	df_dx2 = (f(X[0], X[1]+h/2) - f(X[0], X[1]-h/2))/h

	if n == 0:
		return f(X)

	if n == 1:
		return np.array([-df_dx1, -df_dx2])
	
	if n == 2:
		d2f_dx2 = (f(X[0]+h, X[1]) - 2*f(X[0], X[1])+f(X[0]-h, X[1]))/h**2
		d2f_dy2 = (f(X[0], X[1]+h) - 2*f(X[0], X[1]) + f(X[0], X[1]-h))/h**2
		d2f_dxdy = (f(X[0]+h/2, X[1]+h/2)-f(X[0]+h/2, X[1]-h/2)-f(X[0]-h/2, X[1]+h/2)+f(X[0]-h/2, X[1]-h/2))/h**2

		return np.array([[d2f_dx2, d2f_dxdy], [d2f_dxdy, d2f_dy2]])
    
def phi(f, g, lamb_da, X):
    return f(X[0]+lamb_da*g[0], X[1]+lamb_da*g[1])

def dphi_dlamb_da(lamb_da, n, f, X, d):
    h = 0.01
    
    if n == 1:
        return (phi(f, d, lamb_da+h/2, X) - phi(f, d, lamb_da-h/2, X))/h
    
    if n == 2:
        return (phi(f, d, lamb_da+h, X) - 2*phi(f, d, lamb_da, X) + phi(f, d, lamb_da-h, X))/h**2

def main():
	X0 = np.array([-1.2, 1.0])
	S01 = np.array([1.0, 0.0])
	S02 = np.array([0.0, 1.0])
	#H = dnf_dxn(X0, 2, rosen)

	lamb_da01 = 2.0
	lamb_da02 = 2.0
	lamb_da03 = 2.0
	lamb_da04 = 2.0

	#a = (H[0][0] + np.sqrt(H[0][0]**2 + 4*H[0][1]**2))/2*H[0][1]
	#b = (H[0][1] + np.sqrt(H[0][1]**2 + 4*H[0][0]**2))/2*H[0][0]

	

	for _ in range(100):
		
		for _ in range(20):
			lamb_da1 = lamb_da01 - dphi_dlamb_da(lamb_da01, 1, rosen, X0, S02)/dphi_dlamb_da(lamb_da01, 2, rosen, X0, S02)
			lamb_da01 = lamb_da1
			if dphi_dlamb_da(lamb_da01, 2, rosen, X0, S02) == 0:
				continue
		

		X1 = X0 + lamb_da1*S02

		for _ in range(20):
			lamb_da2 = lamb_da02 - dphi_dlamb_da(lamb_da02, 1, rosen, X1, S01)/dphi_dlamb_da(lamb_da02, 2, rosen, X1, S01)
			lamb_da02 = lamb_da2
			if dphi_dlamb_da(lamb_da01, 2, rosen, X1, S01) == 0:
				continue

		
		

		X02 = X1

		X2 = X02 + lamb_da2*S01

		X03 = X2
		
		for _ in range(20):
			lamb_da3 = lamb_da03 - dphi_dlamb_da(lamb_da03, 1, rosen, X03, S02)/dphi_dlamb_da(lamb_da03, 2, rosen, X03, S02)
			lamb_da03 = lamb_da3
			if dphi_dlamb_da(lamb_da01, 2, rosen, X03, S02) == 0:
				continue

		X3 = X03 + lamb_da3*S02
		
		S = X3 - X1
		X04 = X3
		
		for _ in range(20):
			lamb_da4 = lamb_da04 - dphi_dlamb_da(lamb_da04, 1, rosen, X04, S)/dphi_dlamb_da(lamb_da04, 2, rosen, X04, S)
			lamb_da04 = lamb_da4
			if dphi_dlamb_da(lamb_da01, 2, rosen, X04, S) == 0:
				continue
	

		X4 = X04 + lamb_da4*S
		

		print('X =', X4)

		X0 = X4
		
		
if __name__ == '__main__':
	main()
