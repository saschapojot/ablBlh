from consts import *

#this module is the realization of rk4 algorithm to compute the wavefunction

def x(n):
    return 2 * n * omegaF


def y(n):
    return (2 * n + 1) * omegaF


xmValsAll = [x(m) for m in range(0, N)]
ymValsAll = [y(m) for m in range(0, N)]


def zerothVector(q,aVec, bVec):
    '''

    :param q: time step
    :param aVec: input a vector at time step q, n=0,1,...,N-1
    :param bVec: input b vector at time step q, n=0,1,...,N-1
    :return: K0 vector, M0 vector
    '''
    K0=[]
    M0=[]
    #calculate K0^0:
    K0.append(-1j*v(q*dt)*bVec[0]\
              -1j*u(q*dt)*aVec[0]\
              -1j*xmValsAll[0]*aVec[0]\
              -1j*g*np.abs(aVec[0])**2*bVec[0]
              )
    #calculate K0^n, n=1,2,..., N-1
    for n in range(1,N):
        K0.append(-1j*v(q*dt)*bVec[n]\
                  -1j*w(q*dt)*bVec[n-1]\
                  -1j*u(q*dt)*aVec[n]\
                  -1j*xmValsAll[n]*aVec[n]\
                  -1j*g*np.abs(aVec[n])**2*(bVec[n-1]+bVec[n]))
    #calculate M0^n, n=0,1,...,N-2
    for n in range(0,N-1):
        M0.append(-1j*v(q*dt)*aVec[n]\
                  -1j*w(q*dt)*aVec[n+1]\
                  +1j*u(q*dt)*bVec[n]\
                  -1j*ymValsAll[n]*bVec[n]\
                  -1j*g*np.abs(bVec[n])**2*(aVec[n]+aVec[n+1]))
    #calculate M0^(N-1)
    M0.append(-1j*v(q*dt)*aVec[N-1]\
              +1j*u(q*dt)*bVec[N-1]\
              -1j*ymValsAll[N -1]*bVec[N-1]\
              -1j*g*np.abs(bVec[N-1])**2*aVec[N-1])
    return K0, M0

def firstVecttor(q, aVec, bVec, M0, K0):
    '''

    :param q: time step
    :param aVec: input a vector at time step q, n=0,1,...,N-1
    :param bVec: input b vector at time step q, n=0,1,...,N-1
    :param M0: M0 coef vector
    :param K0: K0 coef vector
    :return: K1, M1
    '''
    K1=[]
    M1=[]
    #calculate K1^0
    K1.append(-1j*v((q+1/2)*dt)*(bVec[0]+1/2*dt*M0[0])\
              -1j*u((q+1/2)*dt)(aVec[0]+1/2*dt*K0[0])\
              -1j*xmValsAll[0]*(aVec[0]+1/2*dt*K0[0])\
              -1j*g*np.abs(aVec[0]+1/2*dt*K0[0])**2*(bVec[0]+1/2*dt*M0[0])

    )
    #calculate K1^n, n=1,2,...,N-1
    for n in range(1,N):
        K1.append(
            -1j*v((q+1/2)*dt)*(bVec[n]+1/2*dt*M0[n])\
            -1j*w((q+1/2)*dt)*(bVec[n-1]+1/2*dt*M0[n-1])\
            -1j*u((q+1/2)*dt)*(aVec[n]+1/2*dt*K0[n])\
            -1j*xmValsAll[n]*(aVec[n]+1/2*dt*K0[n])\
            -1j*g*np.abs(aVec[n]+1/2*dt*K0[n])**2*(bVec[n-1]+1/2*dt*M0[n-1]+bVec[n]+1/2*dt*M0[n])
        )
    #calculate M1^n, n=0,1,...,N-2
    for n in range(0,N-1):
        M1.append(
            -1j*v((q+1/2)*dt)*(aVec[n]+1/2*dt*K0[n])\
            -1j*w((q+1/2)*dt)*(aVec[n+1]+1/2*dt*K0[n+1])\
            +1j*u((q+1/2)*dt)*(bVec[n]+1/2*dt*M0[n])\
            -1j*ymValsAll[n]*(bVec[n]+1/2*dt*M0[n])\
            -1j*g*np.abs(bVec[n]+1/2*dt*M0[n])**2*(aVec[n]+1/2*dt*K0[n]+aVec[n+1]+1/2*dt*K0[n+1])
        )
    #calculate M1^(N-1)
    M1.append(
        -1j*v((q+1/2)*dt)*(aVec[N-1]+1/2*dt*K0[N-1])\
        +1j*u((q+1/2)*dt)*(bVec[N-1]+1/2*dt*M0[N-1])\
        -1j*ymValsAll[N-1]*(bVec[N-1]+1/2*dt*M0[N-1])\
        -1j*g*np.abs(bVec[N-1]+1/2*dt*M0[N-1])**2*(aVec[N-1]+1/2*dt*K0[N-1])
    )

    return K1, M1

def secondVector(q,aVec,bVec,K1,M1):
    '''

    :param q: time step
    :param aVec: input a vector at time step q, n=0,1,...,N-1
    :param bVec: input b vector at time step q, n=0,1,...,N-1
    :param K1: K1 coef vector
    :param M1: M1 coef vector
    :return: K2, M2
    '''
    K2=[]
    M2=[]
    #calculate K2^0
    K2.append(
        -1j*v((q+1/2)*dt)*(bVec[0]+1/2*dt*M1[0])\
        -1j*u((q+1/2)*dt)*(aVec[0]+1/2*dt*K1[0])\
        -1j*xmValsAll[0]*(aVec[0]+1/2*dt*K1[0])\
        -1j*g*np.abs(aVec[0]+1/2*dt*K1[0])**2*(bVec[0]+1/2*dt*M1[0])
    )
    #calculate K2^n, n=1,2,...,N-1
    for n in range(1,N):
        K2.append(
            -1j*v((q+1/2)*dt)*(bVec[n]+1/2*dt*M1[n])\
            -1j*w((q+1/2)*dt)*(bVec[n-1]+1/2*dt*M1[n-1])\
            -1j*u((q+1/2)*dt)*(aVec[n]+1/2*dt*K1[n])\
            -1j*xmValsAll[n]*(aVec[n]+1/2*dt*K1[n])\
            -1j*g*np.abs(aVec[n]+1/2*dt*K1[n])**2*(bVec[n-1]+1/2*dt*M1[n-1]+bVec[n]+1/2*dt*M1[n])
        )
    #calculate M2^n, n=0,1,..., N-2
    for n in range(0,N-1):
        M2.append(
            -1j*v((q+1/2)*dt)*(aVec[n]+1/2*dt*K1[n])\
            -1j*w((q+1/2)*dt)*(aVec[n+1]+1/2*dt*K1[n+1])\
            +1j*u((q+1/2)*dt)*(bVec[n]+1/2*dt*M1[n])\
            -1j*ymValsAll[n]*(bVec[n]+1/2*dt*M1[n])\
            -1j*g*np.abs(bVec[n]+1/2*dt*M1[n])**2*(aVec[n]+1/2*dt*K1[n]+aVec[n+1]+1/2*dt*K1[n+1])
        )
    #calculate M2^(N-1)
    M2.append(
        -1j*v((q+1/2)*dt)*(aVec[N-1]+1/2*dt*K1[N-1])\
        +1j*u((q+1/2)*dt)*(bVec[N-1]+1/2*dt*M1[N-1])\
        -1j*ymValsAll[N-1]*(bVec[N-1]+1/2*dt*M1[N-1])\
        -1j*g*np.abs(bVec[N-1]+1/2*dt*M1[N-1])**2*(aVec[N-1]+1/2*dt*K1[N-1])
    )
    return K2, M2


def thirdVector(q,aVec,bVec,K2, M2):
    '''

    :param q: time step
    :param aVec: input a vector at time step q, n=0,1,...,N-1
    :param bVec: input b vector at time step q, n=0,1,...,N-1
    :param K2: K2 coef vector
    :param M2: M2 coef vector
    :return: K3, M3
    '''
    K3=[]
    M3=[]
    #calculate K3^0
    K3.append(
        -1j*v((q+1)*dt)*(bVec[0]+dt*M2[0])\
        -1j*u((q+1)*dt)*(aVec[0]+dt*K2[0])\
        -1j*xmValsAll[0]*(aVec[0]+dt*K2[0])\
        -1j*g*np.abs(aVec[0]+dt*K2[0])**2*(bVec[0]+dt*M2[0])
    )
    #calculate K3^n, n=1,...,N-1
    for n in range(1,N):
        K3.append(
            -1j*v((q+1)*dt)*(bVec[n]+dt*M2[n])\
            -1j*w((q+1)*dt)*(bVec[n-1]+dt*M2[n-1])\
            -1j*u((q+1)*dt)*(aVec[n]+dt*K2[n])\
            -1j*xmValsAll[n]*(aVec[n]+dt*K2[n])\
            -1j*g*np.abs(aVec[n]+dt*K2[n])**2*(bVec[n-1]+dt*M2[n-1]+bVec[n]+dt*M2[n])
        )
    #calculate M3^n, n=0,1, ..., N-2
    for n in range(0,N-1):
        M3.append(
            -1j*v((q+1)*dt)*(aVec[n]+dt*K2[n])\
            -1j*w((q+1)*dt)*(aVec[n+1]+dt*K2[n+1])\
            +1j*u((q+1)*dt)*(bVec[n]+dt*M2[n])\
            -1j*ymValsAll[n]*(bVec[n]+dt*M2[n])\
            -1j*g*np.abs(bVec[n]+dt*M2[n])**2*(aVec[n]+dt*K2[n]+aVec[n+1]+dt*K2[n+1])
        )
    #calculate M3^(N-1)
    M3.append(
        -1j*v((q+1)*dt)*(aVec[N-1]+dt*K2[N-1])\
        +1j*u((q+1)*dt)*(bVec[N-1]+dt*M2[N-1])\
        -1j*ymValsAll[N-1]*(bVec[N-1]+dt*M2[N-1])\
        -1j*g*np.abs(bVec[N-1]+dt*M2[N-1])**2*(aVec[N-1]+dt*K2[N-1])
    )
    return K3, M3