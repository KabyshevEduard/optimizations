from math import *
import numpy as np

def parab(x,y):
	return x**2 + y**2 - 1.2*x*y

def rosen(x,y):
	return (1-x)**2 + 100 * (y-x**2)**2


def der1(func,m,d,a,h):
        mp = m + np.dot(a+h/2, d)
        mm = m + np.dot(a-h/2, d)
        return (func(mp[0], mp[1]) - func(mm[0], mm[1]))/h

def der2(func,m,d,a,h):
        mp = m + np.dot(a+h, d)
        mm = m + np.dot(a-h, d)
        m = m + np.dot(a, d)
        return (func(mp[0], mp[1]) + func(mm[0], mm[1]) - 2*func(m[0], m[1]))/h**2

def grad(func,m,h):
        x = m[0]
        y = m[1]
        a = (func(x+h/2, y) - func(x-h/2, y))/h
        b = (func(x, y+h/2) - func(x, y-h/2))/h
        return np.array([a, b])

def der1(func,m,d,a,h):
        mp = m + np.dot(a+h/2, d)
        mm = m + np.dot(a-h/2, d)
        return (func(mp[0], mp[1]) - func(mm[0], mm[1]))/h

def der2(func,m,d,a,h):
        mp = m + np.dot(a+h, d)
        mm = m + np.dot(a-h, d)
        m = m + np.dot(a, d)
        return (func(mp[0], mp[1]) + func(mm[0], mm[1]) - 2*func(m[0], m[1]))/h**2

def main():
        startr=np.array([-1.2,1])

        x = startr
        d = -grad(rosen, x, 1e-3)
        xp = np.array([0, 0]) # x на прошлом шаге
        dp = np.array([0, 0]) # d на прошлом шаге

        for i in range(50):
                if(np.dot(grad(rosen, xp, 1e-3),grad(rosen, xp, 1e-3)) != 0):
                        d = -grad(rosen, x, 1e-3) + np.dot(dp, np.dot(grad(rosen, x, 1e-3), grad(rosen, x, 1e-3))/np.dot(grad(rosen, xp, 1e-3), grad(rosen, xp, 1e-3)))
                l = 1.0
                for _ in range(20):
                        if(der2(rosen, x, d, l, 1e-3) == 0):
                                continue
                        l = l - der1(rosen, x, d, l, 1e-3)/der2(rosen, x, d, l, 1e-3)
                xp = x
                dp = d
                x = x + np.dot(l, d)
                print(i, x, 'x функции розенброка')

if __name__ == '__main__':
        main()