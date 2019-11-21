import numpy as np 
from numpy.linalg import inv
import matplotlib.pyplot as plt 

# intializing variables
step_size,alpha,beta,a,b=0.25,1,2,0,1
N=int((b-a)/step_size)
y=np.zeros((N+1,1))
# number of steps for newton ralphson
n=100                                   

# intializng Jacobian_MATRIX (jacobian) with zeros,  Given dirichlet type BVP  Jacobian_MATRIX[0][0]=1 as [g0=y0-alpha] similarly Jacobian_MATRIX[N][N]=1 as gN=yN-beta
Jacobian_MATRIX=np.zeros((N+1,N+1))					 
Jacobian_MATRIX[0][0]=1								 
Jacobian_MATRIX[N][N]=1								 
                   
# intializing y with Linear interpolating formula
for i in range(N+1):
	y[i]=alpha+(beta-alpha)*i/N       


#function calculates value of g at i_th iteration for a given y vector similarly d(y,i),u(y,i),l(y,i)
def g(y,i):                            
	if i==0:
		return y[0]-alpha
	elif i==N:
		return y[N]-beta
	return -(y[i+1]-2*y[i]+y[i-1])+(step_size**2)*(1-((y[i+1]-y[i-1])/2*step_size)**2)

# intializing G(function values) vector with zeros and updating G with intial y vector 
G=np.zeros((N+1,1))                      
for i in range(N+1):					  
	G[i]=g(y,i)

def ldu_calculator(y,i):								 
	return -1+(y[i+1]-y[i-1])/2 ,2, -1-(y[i+1]-y[i-1])*step_size/2

# updating Jacobian_MATRIX with li,di and ui for given i and y
for i in range(1,N):                     
	Jacobian_MATRIX[i][i-1],Jacobian_MATRIX[i][i],Jacobian_MATRIX[i][i+1]=ldu_calculator(y,i)
	

# iterating over n for newton ralphson 
for i in range(n):
    #calculating v by finding inverse of Jacobian_MATRIX & updating y vector by adding v                       
	v=-1*np.dot(inv(Jacobian_MATRIX),G)				
	y=y+v 			
    # calculating new value of Jacobian_MATRIX vector for updated y # calculating new value of G vector for updated y     
	for j in range(1,N):				
		Jacobian_MATRIX[j][j-1],Jacobian_MATRIX[j][j],Jacobian_MATRIX[j][j+1]=ldu_calculator(y,j)
	for j in range(N+1):				
		G[j]=g(y,j)
print(y)
X=[0,0.25,0.5,0.75,1]
Y=[y[i][0] for i in range(len(y))]
plt.plot(X,Y) 
plt.xlabel('x - axis') 
plt.xlabel('x - axis') 
plt.show() 
