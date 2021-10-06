import math
import sys
import numpy as np
from scipy.special import legendre

#alpha is the rotation
alpha=0.02
# 2q is the degree of legendre polynomial to use
q=4
depth=6
gridsize=2**depth+1

thetaGrid=[]
for i in range (gridsize):
    thetaGrid.append(math.pi/(2*(gridsize-1))*i)

#it seems aGrid needs to be much bigger
aGrid=[]
m=2
while m<(2*q+1):
    aGrid.append([])
    j=2
    a=1
    aGrid[int((m-2)/2)].append(a)
    while j<19:
        a=a/(j*(j+2*m+1))
        aGrid[int((m-2)/2)].append(a)
        j+=2
    m+=2

#print(aGrid)

def sigma(bList, x, mu):
    sum=0
    for m in range(q):
        sum+=bList[m]*legendre(2*m+2)(math.cos(mu))*wm(2*m+2, x)
    return (alpha+(1-alpha)*math.sin(x)/x+sum)

def sigmadr(bList, x, mu):
    sum=(1-alpha)*(x*math.cos(x)-math.sin(x))/(x*x)
    for m in range(q):
        sum+=bList[m]*legendre(2*m+2)(math.cos(mu))*wmprime(2*m+2, x)
    return (sum)

def sigmadt(bList, x, mu):
    sum=0
    for m in range(q):
        sum+=bList[m]*np.polyval(legendre(2*m+2).deriv(), math.cos(mu))*wm(2*m+2, x)*(-math.sin(mu))
    return (sum)

def wm(m , x):
    sum=0
    r=0
    while r<19:
        sum+=aGrid[int((m-2)/2)][int(r/2)]*x**r
        r+=2
    return (x**m*sum)

def wmprime(m,x):
    sum1=0
    sum2=0
    r=0
    while r<19:
        sum1+=aGrid[int((m-2)/2)][int(r/2)]*x**r
        r+=2
    r=2
    while r<19:
        sum2+=aGrid[int((m-2)/2)][int(r/2)]*x**(r-1)*r
        r+=2
    return (m*x**(m-1)*sum1+x**m*sum2)

def FindSurface(bList, mu):
    #there may be multiple zeroes which is a problem
    # it is clear that at x=0, the sigma is positive, so we will just increase until we become negative
    x1=0
    x2=0.1
    epsilon=0.0000001
    s=sigma(bList, x2, mu)
    #print(s)
    i=0
    while s>0 and i<100:
        x1+=0.1
        x2+=0.1
        s=sigma(bList, x2, mu)
        i+=1
        print(s)
    #now we have a zero in [x1, x2]
    #implement bisection
    while x2-x1>epsilon and i<1000:
        if sigma(bList,(x1+x2)/2, mu)>0:
            x1=(x1+x2)/2
        else:
            x2=(x1+x2)/2
        i+=1
    return (x1)

def PsiI(bList,x, mu):
    return (-sigma(bList, x, mu)+(1/4)*alpha*x*x*(1-math.cos(mu)*math.cos(mu)))
   
def PsiO(cList, v, x, mu):
    sum=0
    for m in range(q+1):
        sum+=legendre(2*m)(math.cos(mu))*cList[m]/(x**(2*m+1))
    return (v+sum)

def GradI(bList, x, mu):
    sigr=sigmadr(bList, x, mu)
    sigt=sigmadt(bList, x, mu)
    term1=-sigr**2+sigr/2*alpha*x*(1-math.cos(mu)**2)
    term2=math.sin(sigt/x)*math.sin(-sigt/x+alpha*x*x*math.cos(mu)*math.sin(mu))+math.cos(sigt/x)*math.cos(-sigt/x+alpha*x*x*math.cos(mu)*math.sin(mu))
    return(term1*term2)

def GradO(cList, x, mu):
    sum1=0
    sum2=0
    for m in range(q+1):
        sum1+=legendre(2*m)(math.cos(mu))*cList[m]*(2*m+1)/(x**(2*m+2))
        sum2+=np.polyval(legendre(2*m).deriv(), math.cos(mu))*cList[m]*math.sin(mu)/(x**(2*m+2))
    term1=-sigmadr(bList, x, mu)*sum1
    term2=math.sin(sigmadt(bList, x, mu)/x)*math.sin(-sum2)+math.cos(sigmadt(bList, x, mu)/x)*math.cos(-sum2)
    return(term1*term2)

#this is the integration for f. It uses composite simpsons, but we could switch to Romberg and it would probably be better.
def Rombergf(i, surfaceList, bList, cList, v):
    #first we create the necessary matrix
    #this will be designed to be a 6 deep trapezoidal rule Romberg integration
    Rmat=[]
    Rmat.append([math.pi/2*((PsiI(bList, surfaceList[0][1], surfaceList[0][0])-PsiO(cList, v, surfaceList[0][1], surfaceList[0][0]))*legendre(2*i-2)(math.cos(surfaceList[0][0]))*(math.sin(surfaceList[0][0]))+(PsiI(bList, surfaceList[gridsize-1][1], surfaceList[gridsize-1][0])-PsiO(cList, v, surfaceList[gridsize-1][1], surfaceList[gridsize-1][0]))*legendre(2*i-2)(math.cos(surfaceList[gridsize-1][0]))*(math.sin(surfaceList[gridsize-1][0])))])
    for j in range(depth):
        #first append R_{j+1},1
        Rmat.append([])
        sum=0
        step=int((gridsize-1)/2**(j+1))
        for k in range(2**j):
            sum+=(PsiI(bList, surfaceList[(2*k+1)*step][1], surfaceList[(2*k+1)*step][0])-PsiO(cList, v, surfaceList[(2*k+1)*step][1], surfaceList[(2*k+1)*step][0]))*legendre(2*i-2)(math.cos(surfaceList[(2*k+1)*step][0]))*math.sin(surfaceList[(2*k+1)*step][0])
        Rmat[j+1].append(0.5*(Rmat[j][0]+math.pi/(2*2**j)*sum))
        for k in range(j+1):
            Rmat[j+1].append(Rmat[j+1][k]+1/(4**(k+1)-1)*(Rmat[j+1][k]-Rmat[j][k]))
    #print(Rmat)
    return Rmat[depth][depth]

def Rombergfprime(i, surfaceList, bList, cList, v):
    #first we create the necessary matrix
    #this will be designed to be a 6 deep trapezoidal rule Romberg integration
    #(GradI(bList, surfaceList[j][1], surfaceList[j][0])-GradO(cList, surfaceList[j][1], surfaceList[j][0]))
    Rmat=[]
    Rmat.append([math.pi/2*((GradI(bList, surfaceList[0][1], surfaceList[0][0])-GradO(cList,surfaceList[0][1], surfaceList[0][0]))*legendre(2*i-2)(math.cos(surfaceList[0][0]))*(math.sin(surfaceList[0][0]))+(GradI(bList, surfaceList[gridsize-1][1], surfaceList[gridsize-1][0])-GradO(cList, surfaceList[gridsize-1][1], surfaceList[gridsize-1][0]))*legendre(2*i-2)(math.cos(surfaceList[gridsize-1][0]))*(math.sin(surfaceList[gridsize-1][0])))])
    for j in range(depth):
        #first append R_{j+1},1
        Rmat.append([])
        sum=0
        step=int((gridsize-1)/2**(j+1))
        for k in range(2**j):
            sum+=(GradI(bList, surfaceList[(2*k+1)*step][1], surfaceList[(2*k+1)*step][0])-GradO(cList,surfaceList[(2*k+1)*step][1], surfaceList[(2*k+1)*step][0]))*legendre(2*i-2)(math.cos(surfaceList[(2*k+1)*step][0]))*math.sin(surfaceList[(2*k+1)*step][0])
        Rmat[j+1].append(0.5*(Rmat[j][0]+math.pi/(2*2**j)*sum))
        for k in range(j+1):
            Rmat[j+1].append(Rmat[j+1][k]+1/(4**(k+1)-1)*(Rmat[j+1][k]-Rmat[j][k]))
    #print(Rmat)
    return Rmat[depth][depth]


def f(i, surfaceList, bList, cList, v):
    sum=0
    for j in range(0,len(surfaceList)-2,2):
        # it is normally 2h/6 and in burden and faires we get 2h/6=h/3
        length=(1/6)*math.pi/(gridsize-1)
        sum+=length*(PsiI(bList, surfaceList[j][1], surfaceList[j][0])-PsiO(cList, v, surfaceList[j][1], surfaceList[j][0]))*legendre(2*i-2)(math.cos(surfaceList[j][0]))*(math.sin(surfaceList[j][0]))
        sum+=length*4*(PsiI(bList, surfaceList[j+1][1], surfaceList[j+1][0])-PsiO(cList, v, surfaceList[j+1][1], surfaceList[j+1][0]))*legendre(2*i-2)(math.cos(surfaceList[j+1][0]))*(math.sin(surfaceList[j+1][0]))
        sum+=length*(PsiI(bList, surfaceList[j+2][1], surfaceList[j+2][0])-PsiO(cList, v, surfaceList[j+2][1], surfaceList[j+2][0]))*legendre(2*i-2)(math.cos(surfaceList[j+2][0]))*(math.sin(surfaceList[j+2][0]))
        
        #length=(1/6)*((surfaceList[j+2][1])**2+(surfaceList[j][1])**2-2*surfaceList[j][1]*surfaceList[j+2][1]*math.cos(surfaceList[j+2][0]-surfaceList[j][0]))**(1/2)
        #sum+=length*(PsiI(bList, surfaceList[j][1], surfaceList[j][0])-PsiO(cList, v, surfaceList[j][1], surfaceList[j][0]))*legendre(2*i-2)(math.cos(surfaceList[j][0]))
        #sum+=length*4*(PsiI(bList, surfaceList[j+1][1], surfaceList[j+1][0])-PsiO(cList, v, surfaceList[j+1][1], surfaceList[j+1][0]))*legendre(2*i-2)(math.cos(surfaceList[j+1][0]))
        #sum+=length*(PsiI(bList, surfaceList[j+2][1], surfaceList[j+2][0])-PsiO(cList, v, surfaceList[j+2][1], surfaceList[j+2][0]))*legendre(2*i-2)(math.cos(surfaceList[j+2][0]))
    return 2*sum

#this is the integration for f. It uses composite simpsons, but we could switch to Romberg and it would probably be better.
def fprime(i, surfaceList, bList, cList, v):
    sum=0
    #I have changed the formula here. It seems using the real length accounted for the sin at the end of the integrand
    for j in range(0,len(surfaceList)-2, 2):

        length=(1/6)*math.pi/(gridsize-1)
        sum+=length*(GradI(bList, surfaceList[j][1], surfaceList[j][0])-GradO(cList, surfaceList[j][1], surfaceList[j][0]))*legendre(2*i-4)(math.cos(surfaceList[j][0]))*(math.sin(surfaceList[j][0]))
        sum+=length*4*(GradI(bList, surfaceList[j+1][1], surfaceList[j+1][0])-GradO(cList, surfaceList[j+1][1], surfaceList[j+1][0]))*legendre(2*i-4)(math.cos(surfaceList[j+1][0]))*(math.sin(surfaceList[j+1][0]))
        sum+=length*(GradI(bList, surfaceList[j+2][1], surfaceList[j+2][0])-GradO(cList, surfaceList[j+2][1], surfaceList[j+2][0]))*legendre(2*i-4)(math.cos(surfaceList[j+2][0]))*(math.sin(surfaceList[j+2][0]))
        
        #length=(1/6)*((surfaceList[j+2][1])**2+(surfaceList[j][1])**2-2*surfaceList[j][1]*surfaceList[j+2][1]*math.cos(surfaceList[j+2][0]-surfaceList[j][0]))**(1/2) 
        #sum+=length*(GradI(bList, surfaceList[j][1], surfaceList[j][0])-GradO(cList, surfaceList[j][1], surfaceList[j][0]))*legendre(2*i-4)(math.cos(surfaceList[j][0]))
        #sum+=length*4*(GradI(bList, surfaceList[j+1][1], surfaceList[j+1][0])-GradO(cList, surfaceList[j+1][1], surfaceList[j+1][0]))*legendre(2*i-4)(math.cos(surfaceList[j+1][0]))
        #sum+=length*(GradI(bList, surfaceList[j+2][1], surfaceList[j+2][0])-GradO(cList, surfaceList[j+2][1], surfaceList[j+2][0]))*legendre(2*i-4)(math.cos(surfaceList[j+2][0]))
    return 2*sum

# now we can begin the iteration
# we will start with b equal to zero and c,v =1.
#if q is 7, we have 7 b's and 8 c's

#bList cannot start at zero because of how this guy did the things. But you must start very very small
bList=[-0.013, -0.00001, -0.0000001, -0.00000001]
v=1
cList=[]
for i in range(q):
    #bList.append(0.0001)
    cList.append(1)
cList.append(1)

surfaceList=[]

print(FindSurface( bList, math.pi/2) )
#for mu in thetaGrid:
#    surfaceList.append([mu, FindSurface(bList, mu)])
#print(surfaceList)
#print(Rombergf(i+1,surfaceList, bList, cList, v))

maxdiff=1
iter=1
print (v)
print(bList)
print(cList)


while (maxdiff>0.01 and iter<1):
    
    #with our guess, the first thing we do is compute the surface estimate
    #for Romberg integration, we do not need to use the same thetaGrid. 
    #it appears that at the level they want, Romberg would only need to use 33 points? 

    surfaceList=[]
    for mu in thetaGrid:
        surfaceList.append([mu, FindSurface(bList, mu)])

    fname="iteration{}surface.txt".format(iter)  
    with open (fname, "w") as file:
        for i in range(len(surfaceList)):
            file.write(str(surfaceList[i][0])+" "+str(surfaceList[i][1]))
            file.write("\n")
    file.close()

    #print(surfaceList)

    #the (v, bList, cList)=x in Williams' paper. In the iteration, the change in x is given by solving G(dx)=-f
    # we compute f

    fList=[]
    for i in range(1, q+2):
        fList.append(-Rombergf(i, surfaceList, bList, cList, v))


    for i in range(2, q+3):
        fList.append(-Rombergfprime(i, surfaceList, bList, cList, v))

    #print (fList)
    #now we compute G
    s=0.001
    Jacobian=[]

    for i in range(q+1):
        Jacobian.append([])
        #note we follow williams here on the derivative estimates, but it is unclear if these are quite right. It could be used to avoid subtraction of small numbers.
        # it seems accuracy issues are due to derivative estimates 
        h=v*s
        print(Rombergf(i+1,surfaceList, bList, cList, v+h))
        Jacobian[i].append(Rombergf(i+1,surfaceList, bList, cList, v+h)/(h))
        for j in range(q):
            hplusList=[]
            hminusList=[]
            for k in range(len(bList)):
                if k==j:
                    h=bList[j]*s
                    hplusList.append(bList[j]+h)
                    hminusList.append(bList[j]-h)
                else:
                    hplusList.append(bList[j])
                    hminusList.append(bList[j])
            Jacobian[i].append((Rombergf(i+1, surfaceList, hplusList, cList, v))/(h))
        for j in range (q+1):
            hplusList=[]
            hminusList=[]
            for k in range(len(cList)):
                if k==j:
                    h=cList[j]*s
                    hplusList.append(cList[j]+h)
                    hminusList.append(cList[j]-h)
                else:
                    hplusList.append(cList[j])
                    hminusList.append(cList[j])
            Jacobian[i].append((Rombergf(i+1, surfaceList, bList, hplusList, v))/(h))
        #print(Jacobian[i])
    for i in range(2,q+3):
        Jacobian.append([])
        h=v*s
        Jacobian[i+q-1].append((Rombergfprime(i,surfaceList, bList, cList, v+h))/(h))

        for j in range(q):
            hplusList=[]
            hminusList=[]
            for k in range(len(bList)):
                if k==j:
                    h=bList[j]*s
                    hplusList.append(bList[j]+h)
                    hminusList.append(bList[j]-h)
                else:
                    hplusList.append(bList[j])
                    hminusList.append(bList[j])
            Jacobian[i+q-1].append((Rombergfprime(i, surfaceList, hplusList, cList, v))/(h))
        for j in range (q+1):
            hplusList=[]
            hminusList=[]
            for k in range(len(cList)):
                if k==j:
                    h=cList[j]*s
                    hplusList.append(cList[j]+h)
                    hminusList.append(cList[j]-h)
                else:
                    hplusList.append(cList[j])
                    hminusList.append(cList[j])
            Jacobian[i+q-1].append((Rombergfprime(i, surfaceList, bList, hplusList, v))/(h))
        #print(Jacobian[i+q-1])
    print(Jacobian)
    soln=np.linalg.solve(Jacobian, fList)
    maxdiff=0
    for i in range(2*q+2):
        if abs(soln[i])>maxdiff:
            maxdiff=abs(soln[i])
    
    if maxdiff<0.01:
        print('Success!')
    v+=soln[0]
    for i in range(len(bList)):
        bList[i]+=soln[1+i]
    for i in range(len(cList)):
        cList[i]+=soln[1+i+len(bList)]
    print (v)
    print(bList)
    print(cList)
    iter+=1