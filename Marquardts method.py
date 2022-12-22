import numpy as np

def f(x1, x2):
	return x1**2 + x2**2 - 1.2*x1*x2

def rosen(x1, x2):
	return (1-x1)**2 + 100 * (x2 - x1**2)**2
 
def dnf_dxn(X, n, f):
	h = 0.0001
	df_dx1 = (f(X[0]+h/2, X[1]) - f(X[0]-h/2, X[1]))/h
	df_dx2 = (f(X[0], X[1]+h/2) - f(X[0], X[1]-h/2))/h

	if n == 0:
		return f(X)

	if n == 1:
		return np.array([-df_dx1, -df_dx2])
	
	else:
		d2f_dx2 = (f(X[0]+h, X[1]) - 2*f(X[0], X[1])+f(X[0]-h, X[1]))/h**2
		d2f_dy2 = (f(X[0], X[1]+h) - 2*f(X[0], X[1]) + f(X[0], X[1]-h))/h**2
		d2f_dxdy = (f(X[0]+h/2, X[1]+h/2)-f(X[0]+h/2, X[1]-h/2)-f(X[0]-h/2, X[1]+h/2)+f(X[0]-h/2, X[1]-h/2))/h**2

		return np.array([[d2f_dx2, d2f_dxdy], [d2f_dxdy, d2f_dy2]])

def main():
	X = np.array([4.0, 1.0])
	E = np.array([[1.0, 0], [0, 1.0]])
	alpha = 4.0
	betta = 0.5

	for _ in range(100):
		H = dnf_dxn(X, 2, rosen)
		g = dnf_dxn(X, 1, rosen)
		S = np.dot(np.linalg.inv(np.add(H, np.multiply(alpha, E))), g)

		X = X + S
		print('X =', X)

		alpha = alpha * betta

if __name__ == '__main__':
	main()
