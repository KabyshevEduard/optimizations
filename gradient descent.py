import numpy as np
from scipy.special import logsumexp

def f(x1, x2):
	return x1**2 + x2**2 - 1.2*x1*x2


def phi(x1, x2, alpha, grad):
	
	a1 = x1 - alpha * grad[0]
	a2 = x2 - alpha * grad[1]

	return f(a1, a2)


def grad(h, X):

	a = (f(X[0] + h/2, X[1]) - f(X[0] - h/2, X[1]))/h
	b = (f(X[0], X[1] + h/2) - f(X[0], X[1] - h/2))/h
	g1 = logsumexp(a/np.sqrt(a**2 + b**2))
	g2 = logsumexp(b/np.sqrt(a**2 + b**2))

	return np.array([g1, g2]).astype(np.float32)


def Phi_alpha(h, X, alpha):

	
	x = X[0]
	y = X[1]
	g = grad(h, X)
	phi_alpha = (phi(x, y, alpha + h/2, g) - phi(x, y, alpha - h/2, g))/h

	return phi_alpha


def Phi_2alpha(h, X, alpha):

	x = X[0]
	y = X[1]
	g = grad(h, X)
	phi_2alpha = (phi(x, y, alpha + h, g) - 2*phi(x, y, alpha, g) + phi(x, y, alpha - h, g))/h**2
	
	return phi_2alpha


alpha0 = 0.0
X = np.array([4.0, 1.0]).astype(np.float32)
h = 0.001



for i in range(61):
	j = grad(h, X)
	print('j = ', j)

	alpha = alpha0 - Phi_alpha(h, X, alpha0)/Phi_2alpha(h, X, alpha0)
	X[0] = X[0] - alpha*j[0]
	X[1] = X[1] - alpha*j[1]
	print('X = ', X)
	print('lambda = ', alpha)

	alpha0 = alpha
	

	