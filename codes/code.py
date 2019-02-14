'''
Question ->
The normal to the curve, x^2 + 2xy - 3y^2 = 0, at (1, 1):
(1) does not meet the curve again.
(2)  meets  the  curve  again  in  the  second quadrant.
(3)  meets  the  curve  again  in  the  third quadrant.
(4)  meets  the  curve  again  in  the  fourth quadrant.
'''

import numpy as np
import matplotlib.pyplot as plt

#x^2 + 2xy - 3y^2 = 0 passes through origin

#returns a point such that x coordinates of P_new and P differ by del_x and the line joining P to P_new has the direction vector m
def point(P, m, del_x):		
	P_new = P + m / m[0] * del_x  
	return P_new
	
#draws a line segment joining points A and B	
def line_seg(A, B, name, c):	
	size = int(abs((B[0]-A[0])*100))
	l = np.linspace(0, 1, size)
	o = np.ones(size)
	x_i = np.matmul(np.reshape(A, (2,1)), np.reshape(o, (1, size))) + np.matmul(np.reshape((B-A), (2,1)), np.reshape(l, (1, size))) #A + l(B-A)
	M = B - A
	plt.plot(x_i[0,:], x_i[1,:], label = (name + ': y = '+str(round(M[1]/M[0],2))+' x + '+str(c)))
	plt.plot(A[0], A[1], 'o')
	plt.plot(B[0], B[1], 'o')	
	
#returns direction vectors of the lines in the combined line equation ax^2 + bx + c = 0 in the form of a matrix (m1 m2), where m = (x y)^T	
def separate_lines(a, b, c): 	
	M = np.ones((2, 2))
	M[0][0] = (-b + (b**2 - 4*a*c)**0.5)/(2*a)
	M[0][1] = (-b - (b**2 - 4*a*c)**0.5)/(2*a)
	return M
	
#given a pair of lines, if line has normal with d.v. m and passes through P, finding X, such that X satisfies (X-P) x m = 0 for both lines
#X [m1 m2] = [P1xm1 P2xm2] => X M = K
def find_intersection(m1, m2, P1, P2): 
	M = np.transpose(np.vstack((m1, m2)))
	M_inv = np.linalg.inv(M)
	K = np.array([np.dot(P1, m1), np.dot(P2, m2)])
	X = np.matmul(K, M_inv)
	return X

	
#plotting origin and (1,1)
O = np.array([0, 0]); P = np.array([1, 1])
plt.plot(0, 0, 'o')
plt.plot(1, 1, 'o')
plt.annotate('O = '+str(O), xy=(O[0], O[1]), xytext=(O[0], O[1]-1))
plt.annotate('P = '+str(P), xy=(P[0], P[1]), xytext=(P[0]+0.3, P[1]))

#setting parameters of the combined line equation
a = 1; b = 2; c = -3
M = separate_lines(a, b, c)

#drawing the lines
for i in range(0, 2):
	A = point(O, np.transpose(M)[i], 5)
	B = point(O, np.transpose(M)[i], -5)
	line_seg(A, B, 'Line '+str(i+1), 0)
	
#N stores the d.v.s of the normals of the lines
N = np.dot(np.array([[0, -1], [1, 0]]), M)
	
#n stores index of line on which (1,1) lies and n_d the other
n = int(np.dot(N[1], P) == 0)
n_d = int(not bool(n))

#m stores its d.v. 
m = np.transpose(N)[n]
A = point(P, m, 5)
B = point(P, m, -5)
c = float(np.dot(np.transpose(M)[n], P))/m[1]
line_seg(A, B, 'Normal to line '+str(n+1), c)

#checking if the lines intersect before finding out and plotting the point of intersection
if float(M[1][n])/M[0][n] != float(N[1][n_d])/N[0][n_d]: #checking if normals of the two lines have same d.v.
	I = find_intersection(np.transpose(M)[n], np.transpose(N)[n_d], P, O)
	plt.plot(I[0], I[1], 'o')
	plt.annotate('I = '+str(I), xy=(I[0], I[1]), xytext=(I[0]+0.1, I[1]+0.1))
	
plt.axis('equal')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc = 'upper right')
plt.grid() 
plt.show()

'''
Answer from graph -> (4) meets  the  curve  again  in  the  fourth quadrant.
'''
