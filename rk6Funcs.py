from consts import *

#this module is the realization of rk6 algorithm to compute the wavefunction

def x(n):
    return 2 * n * omegaF


def y(n):
    return (2 * n + 1) * omegaF


xmValsAll = [x(m) for m in range(0, N)]
ymValsAll = [y(m) for m in range(0, N)]


def firstVector(q,h,g,aVec,bVec):
    '''

    :param q: time step
    :param h: time duration
    :param g: nonlinearity
    :param aVec: input wvfcnt at even sites
    :param bVec: input wvfcnt at odd sites
    :return: K1, M1
    '''
    K1=[]
    M1=[]
    K1.append(
        -1j*v(q*h)*bVec[0]-1j*u(q*h)*aVec[0]-1j*xmValsAll[0]*aVec[0]-1j*g*np.abs(aVec[0])**2*bVec[0]
    )
    #n=1,...,N-1
    for n in range(1,N):
        K1.append(
            -1j*v(q*h)*bVec[n]-1j*w(q*h)*bVec[n-1]-1j*u(q*h)*aVec[n]-1j*xmValsAll[n]*aVec[n]\
        -1j*g*np.abs(aVec[n])**2*(bVec[n-1]+bVec[n])
        )

    #n=0,...,N-2
    for n in range(0,N-1):
        M1.append(
            -1j*v(q*h)*aVec[n]-1j*w(q*h)*aVec[n+1]+1j*u(q*h)*bVec[n]-1j*ymValsAll[n]*bVec[n]\
            -1j*g*np.abs(bVec[n])**2*(aVec[n]+aVec[n+1])
        )
    M1.append(
        -1j*v(q*h)*aVec[N-1]+1j*u(q*h)*bVec[N-1]-1j*ymValsAll[N-1]*bVec[N-1]\
        -1j*g*np.abs(bVec[N-1])**2*aVec[N-1]
    )
    return K1, M1



def secondVector(q,h,g,aVec,bVec,K1,M1):
    '''

    :param q: time step
    :param h: time duration
    :param g: nonlinearity
    :param aVec: input wvfcnt at even sites
    :param bVec:  input wvfcnt at odd sites
    :param K1: K1 coef vector
    :param M1: M1 coef vector
    :return: K2, M2
    '''
    K2=[]
    M2=[]
    K2.append(
        -1j*v(q*h+h)*(bVec[0]+h*M1[0])-1j*u(q*h+h)*(aVec[0]+h*K1[0])\
        -1j*xmValsAll[0]*(aVec[0]+h*K1[0])\
        -1j*g*np.abs(aVec[0]+h*K1[0])**2*(bVec[0]+h*M1[0])

    )
    #n=1,...,N-1
    for n in range(1,N):
        K2.append(
            -1j*v(q*h+h)*(bVec[n]+h*M1[n])-1j*w(q*h+h)*(bVec[n-1]+h*M1[n-1])\
            -1j*u(q*h+h)*(aVec[n]+h*K1[n])-1j*xmValsAll[n]*(aVec[n]+h*K1[n])\
            -1j*g*np.abs(aVec[n]+h*K1[n])**2*(bVec[n-1]+h*K1[n-1]+bVec[n]+h*K1[n])
        )

    #n=0,...,N-2
    for n in range(0,N-1):
        M2.append(
            -1j*v(q*h+h)*(aVec[n]+h*K1[n])-1j*w(q*h+h)*(aVec[n+1]+h*K1[n+1])\
            +1j*u(q*h+h)*(bVec[n]+h*M1[n])-1j*ymValsAll[n]*(bVec[n]+h*M1[n])\
            -1j*g*np.abs(bVec[n]+h*M1[n])**2*(aVec[n]+h*K1[n]+aVec[n+1]+h*K1[n+1])
        )
    M2.append(
        -1j*v(q*h+h)*(aVec[N-1]+h*K1[N-1])+1j*u(q*h+h)*(bVec[N-1]+h*M1[N-1])\
        -1j*ymValsAll[N-1]*(bVec[N-1]+h*M1[N-1])\
        -1j*g*np.abs(bVec[N-1]+h*M1[N-1])**2*(aVec[N-1]+h*K1[N-1])
    )

    return K2, M2


def thirdVector(q,h,g,aVec,bVec,K1,M1,K2,M2):
    '''

    :param q: time step
    :param h: time duration
    :param g: nonlinearity
    :param aVec: input wvfcnt at even sites
    :param bVec: input wvfcnt at odd sites
    :param K1: K1 coef vector
    :param M1: M1 coef vector
    :param K2: K2 coef vector
    :param M2: M2 coef vector
    :return: K3, M3
    '''
    K3=[]
    M3=[]
    timeQuad=q*h+1/2*h
    K3.append(
        -1j*v(timeQuad)*(bVec[0]+3/8*h*M1[0]+1/8*h*M2[0])\
        -1j*u(timeQuad)*(aVec[0]+3/8*h*K1[0]+1/8*h*K2[0])\
        -1j*xmValsAll[0]*(aVec[0]+3/8*h*K1[0]+1/8*h*K2[0])
    )
    #n=1,...,N-1
    for n in range(1,N):
        K3.append(
            -1j*v(timeQuad)*(bVec[n]+3/8*h*M1[n]+1/8*h*M2[n])\
            -1j*w(timeQuad)*(bVec[n-1]+3/8*h*M1[n-1]+1/8*h*M2[n-1])\
            -1j*u(timeQuad)*(aVec[n]+3/8*h*K1[n]+1/8*h*K2[n])\
            -1j*xmValsAll[n]*(aVec[n]+3/8*h*K1[n]+1/8*h*K2[n])\
            -1j*g*np.abs(aVec[n]+3/8*h*K1[n]+1/8*h*K2[n])**2\
            *(bVec[n-1]+3/8*h*M1[n-1]+1/8*h*M2[n-1]+bVec[n]+3/8*h*M1[n]+1/8*h*M2[n])
        )
    #n=0,...,N-2
    for n in range(0,N-1):
        M3.append(
            -1j*v(timeQuad)*(aVec[n]+3/8*h*K1[n]+1/8*h*K2[n])\
            -1j*w(timeQuad)*(aVec[n+1]+3/8*h*K1[n+1]+1/8*h*K2[n+1])\
            +1j*u(timeQuad)*(bVec[n]+3/8*h*M1[n]+1/8*h*M2[n])\
            -1j*ymValsAll[n]*(bVec[n]+3/8*h*M1[n]+1/8*h*M2[n])\
            -1j*g*np.abs(bVec[n]+3/8*h*M1[n]+1/8*h*M2[n])**2\
            *(aVec[n]+3/8*h*K1[n]+1/8*h*K2[n]+aVec[n+1]+3/8*h*K1[n+1]+1/8*h*K2[n+1])
        )
    M3.append(
        -1j*v(timeQuad)*(aVec[N-1]+3/8*h*K1[N-1]+1/8*h*K2[N-1])\
        +1j*u(timeQuad)*(bVec[N-1]+3/8*h*M1[N-1]+1/8*h*M2[N-1])\
        -1j*ymValsAll[N-1]*(bVec[N-1]+3/8*h*M1[N-1]+1/8*h*M2[N-1])\
        -1j*g*np.abs(bVec[N-1]+3/8*h*M1[N-1]+1/8*h*M2[N-1])**2\
        *(aVec[N-1]+3/8*h*K1[N-1]+1/8*h*K2[N-1])
    )

    return K3, M3



def fourthVector(q,h,g,aVec,bVec,K1,M1,K2,M2,K3,M3):
    '''

    :param q: time step
    :param h: time duration
    :param g: nonlinearity
    :param aVec: input wvfcnt at even sites
    :param bVec: input wvfcnt at odd sites
    :param K1: K1 coef vector
    :param M1: M1 coef vector
    :param K2: K2 coef vector
    :param M2: M2 coef vector
    :param K3: K3 coef vector
    :param M3: M3 coef vector
    :return: K4, M4
    '''
    K4=[]
    M4=[]
    timeQuad=q*h+2/3*h
    K4.append(
        -1j*v(timeQuad)*(bVec[0]+8/27*h*M1[0]+2/27*h*M2[0]+8/27*h*M3[0])\
        -1j*u(timeQuad)*(aVec[0]+8/27*h*K1[0]+2/27*h*K2[0]+8/27*h*K3[0])\
        -1j*xmValsAll[0]*(aVec[0]+8/27*h*K1[0]+2/27*h*K2[0]+8/27*h*K3[0])\
        -1j*g*np.abs(aVec[0]+8/27*h*K1[0]+2/27*h*K2[0]+8/27*h*K3[0])**2\
        *(bVec[0]+8/27*h*M1[0]+2/27*h*M2[0]+8/27*h*M3[0])
    )
    #n=1,...,N-1
    for n in range(1,N):
        K4.append(
            -1j*v(timeQuad)*(bVec[n]+8/27*h*M1[n]+2/27*h*M2[n]+8/27*h*M3[n])\
            -1j*w(timeQuad)*(bVec[n-1]+8/27*h*M1[n-1]+2/27*h*M2[n-1]+8/27*h*M3[n-1])\
            -1j*u(timeQuad)*(aVec[n]+8/27*h*K1[n]+2/27*h*K2[n]+8/27*h*K3[n])\
            -1j*xmValsAll[n]*(aVec[n]+8/27*h*K1[n]+2/27*h*K2[n]+8/27*h*K3[n])\
            -1j*g*np.abs(aVec[n]+8/27*h*K1[n]+2/27*h*K2[n]+8/27*h*K3[n])**2\
            *(bVec[n-1]+8/27*h*M1[n-1]+2/27*h*M2[n-1]+8/27*h*M3[n-1]
              +bVec[n]+8/27*h*M1[n]+2/27*h*M2[n]+8/27*h*M3[n])
        )
    #n=0,...,N-2
    for n in range(0,N-1):
        M4.append(
            -1j*v(timeQuad)*(aVec[n]+8/27*h*K1[n]+2/27*h*K2[n]+8/27*h*K3[n])\
            -1j*w(timeQuad)*(aVec[n+1]+8/27*h*K1[n+1]+2/27*h*K2[n+1]+8/27*h*K3[n+1])\
            +1j*u(timeQuad)*(bVec[n]+8/27*h*M1[n]+2/27*h*M2[n]+8/27*h*M3[n])\
            -1j*ymValsAll[n]*(bVec[n]+8/27*h*M1[n]+2/27*h*M2[n]+8/27*h*M3[n])\
            -1j*g*np.abs(bVec[n]+8/27*h*M1[n]+2/27*h*M2[n]+8/27*h*M3[n])**2\
            *(aVec[n]+8/27*h*K1[n]+2/27*h*K2[n]+8/27*h*K3[n]
              +aVec[n+1]+8/27*h*K1[n+1]+2/27*h*K2[n+1]+8/27*h*K3[n+1])
        )
    M4.append(
        -1j*v(timeQuad)*(aVec[N-1]+8/27*h*K1[N-1]+2/27*h*K2[N-1]+8/27*h*K3[N-1])\
        +1j*u(timeQuad)*(bVec[N-1]+8/27*h*M1[N-1]+2/27*h*M2[N-1]+8/27*h*M3[N-1])\
        -1j*ymValsAll[N-1]*(bVec[N-1]+8/27*h*M1[N-1]+2/27*h*M2[N-1]+8/27*h*M3[N-1])\
        -1j*g*np.abs(bVec[N-1]+8/27*h*M1[N-1]+2/27*h*M2[N-1]+8/27*h*M3[N-1])**2\
        *(aVec[N-1]+8/27*h*K1[N-1]+2/27*h*K2[N-1]+8/27*h*K3[N-1])
    )

    return K4, M4