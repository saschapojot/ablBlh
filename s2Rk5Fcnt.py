from consts import *


def x(n):
    return 2 * n * omegaF


def y(n):
    return (2 * n + 1) * omegaF


xmValsAll = [x(m) for m in range(0, N)]
ymValsAll = [y(m) for m in range(0, N)]
def phi1A2(deltaTau, uQ, alphaVec):
    '''

    :param deltaTau: time step
    :param uQ: u value at time (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :return:
    '''
    for m in range(0, N):
        alphaVec[m] *= np.exp(-1j * deltaTau * (xmValsAll[m] + uQ))


def phi1B2(deltaTau, uQ, betaVec):
    '''

    :param deltaTau: time step
    :param uQ: u value at time step (q+1/2)dt
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    for m in range(0, N):
        betaVec[m] *= np.exp(-1j * deltaTau * (ymValsAll[m] - uQ))


def phi1C1(deltaTau, vQ, alphaVec, betaVec):
    '''

    :param deltaTau: time step
    :param vQ: v value at time step (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    for m in range(0, N):
        alphaVec[m] -= 1j * deltaTau * vQ * betaVec[m]


def phi1D1(deltaTau, vQ, alphaVec, betaVec):
    '''

    :param deltaTau: time step
    :param vQ: v value at time step (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    for m in range(0, N):
        betaVec[m] -= 1j * deltaTau * vQ * alphaVec[m]


def phi21(deltaTau, wQ, alphaVec, betaVec):
    '''

    :param deltaTau: time step
    :param wQ: w value at (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    for j in range(1, N):
        alphaVec[j] -= 1j * deltaTau * wQ * betaVec[j - 1]


def phi22(deltaTau, wQ, alphaVec, betaVec):
    '''

    :param deltaTau: time step
    :param wQ: w value at time step (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    for j in range(0, N - 1):
        betaVec[j] -= 1j * deltaTau * wQ * alphaVec[j + 1]


def zeta1(uQ, vQ, alphaVec, betaVec):
    '''
    composite mapping zeta1
    :param uQ: u value at time (q+1/2)dt
    :param vQ: v value at time (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    phi1B2(dt / 2, uQ, betaVec)
    phi1A2(dt / 2, uQ, alphaVec)
    phi1C1(dt / 2, vQ, alphaVec, betaVec)
    phi1D1(dt, vQ, alphaVec, betaVec)
    phi1C1(dt / 2, vQ, alphaVec, betaVec)
    phi1B2(dt / 2, uQ, betaVec)
    phi1A2(dt / 2, uQ, alphaVec)


def zeta2(wQ, alphaVec, betaVec):
    '''
    composite mapping zeta2
    :param wQ: w value at time step (q+1/2)dt
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''
    phi21(dt / 4, wQ, alphaVec, betaVec)
    phi22(dt / 2, wQ, alphaVec, betaVec)
    phi21(dt / 4, wQ, alphaVec, betaVec)


def expmidtF(q, alphaVec, betaVec):
    '''
    linear step in S2
    :param q: time step q
    :param alphaVec: input wvfcnt at even sites
    :param betaVec: input wvfcnt at odd sites
    :return:
    '''

    t = (q + 1 / 2) * dt
    uq = u(t)
    vq = v(t)
    wq = w(t)
    zeta2(wq, alphaVec, betaVec)
    zeta1(uq, vq, alphaVec, betaVec)
    zeta2(wq, alphaVec, betaVec)

def meanX(psiQ):
    '''

    :param psiQ: wvfcnt at time step q
    :return: mean position at time q
    '''
    xOut=0
    norm2=0
    for j in range(0,len(psiQ)):
        xOut+=j*np.abs(psiQ[j])**2
        norm2+=np.abs(psiQ[j])**2
    return xOut/norm2


def firstVecRk5(firstVec,secondVec,g):
    '''

    :param firstVec: wvfcnt at even sites
    :param secondVec: wvfcnt at odd sites
    :param g: nonlinearity
    :param h: time step
    :return: K1, M1
    '''
    K1=[]
    M1=[]
    K1.append(
        -1j*g*np.abs(firstVec[0])**2*secondVec[0]
    )
    for n in range(1,N):
        #n=1,...,N-1
        K1.append(
            -1j*g*np.abs(firstVec[n])**2*(secondVec[n-1]+secondVec[n])
        )
    for n in range(0,N-1):
        #n=0,...,N-2
        M1.append(
            -1j*g*np.abs(secondVec[n])**2*(firstVec[n]+firstVec[n+1])
        )
    M1.append(
        -1j*np.abs(secondVec[N-1])**2*firstVec[N-1]
    )
    return K1, M1


def secondVecRk5(firstVec,secondVec,g,h,K1,M1):
    '''

    :param firstVec: wvfcnt at even sites
    :param secondVec: wvfcnt at odd sites
    :param g: nonlinearity
    :param h: time step
    :param K1: K1 coef vector
    :param M1: M1 coef vector
    :return: K2, M2
    '''
    K2=[]
    M2=[]
    K2.append(
        -1j*g*np.abs(firstVec[0]+1/4*h*K1[0])**2*(secondVec[0]+1/4*h*M1[0])
    )
    #n=1,...,N-1
    for n in range(1,N):
        K2.append(
            -1j*g*np.abs(firstVec[n]+1/4*h*K1[n])**2\
            *(secondVec[n-1]+1/4*h*M1[n-1]+secondVec[n]+1/4*h*M1[n])
        )
    #n=0,...,N-2
    for n in range(0,N-1):
        M2.append(
            -1j*g*np.abs(secondVec[n]+1/4*h*M1[n])**2\
            *(firstVec[n]+1/4*h*K1[n]+firstVec[n+1]+1/4*h*K1[n+1])
        )
    M2.append(
        -1j*g*np.abs(secondVec[N-1]+1/4*h*M1[N-1])**2\
        *(firstVec[N-1]+1/4*h*K1[N-1])
    )
    return K2, M2

def thirdVecRk5(firstVec,secondVec,g,h,K1,M1,K2,M2):
    '''

    :param firstVec: wvfcnt at even sites
    :param secondVec: wvfcnt at odd sites
    :param g: nonlinearity
    :param h: time step
    :param K1: K1 coef vector
    :param M1: M1 coef vector
    :param K2: K2 coef vector
    :param M2: M2 coef vector
    :return: K3, M3
    '''
    K3=[]
    M3=[]
    K3.append(
        -1j*g*np.abs(firstVec[0]+3/32*h*K1[0]+9/32*h*K2[0])**2\
        *(secondVec[0]+3/32*h*M1[0]+9/32*h*M2[0])
    )
    #n=1,...,N-1
    for n in range(1,N):
        K3.append(
            -1j*g*np.abs(firstVec[n]+3/32*h*K1[n]+9/32*h*K2[n])**2\
            *(secondVec[n-1]+3/32*h*M1[n-1]+9/32*h*M2[n-1]
              +secondVec[n]+3/32*h*M1[n]+9/32*h*M2[n])
        )
    #n=0,...,N-2
    for n in range(0,N-1):
        M3.append(
            -1j*g*np.abs(secondVec[n]+3/32*h*M1[n]+9/32*h*M2[n])**2\
            *(firstVec[n]+3/32*h*K1[n]+9/32*h*K2[n]
              +firstVec[n+1]+3/32*h*K1[n+1]+9/32*h*K2[n+1])
        )
    M3.append(
        -1j*g*np.abs(secondVec[N-1]+3/32*h*M1[N-1]+9/32*h*M2[N-1])**2\
        *(firstVec[N-1]+3/32*h*K1[N-1]+9/32*h*K2[N-1])
    )

    return K3, M3


def fourthVecRk5(firstVec,secondVec,g,h,K1,M1,K2,M2,K3,M3):
    '''

    :param firstVec: wvfcnt at even sites
    :param secondVec: wvfcnt at odd sites
    :param g: nonlinearity
    :param h: time step
    :param K1: K1 coef vector
    :param M1: M1 coef vector
    :param K2: K2 coef vetcor
    :param M2: M2 coef vector
    :param K3: K3 coef vector
    :param M3: M3 coef vecror
    :return: K4, M4
    '''
    K4=[]
    M4=[]
    K4.append(
        -1j*g*np.abs(firstVec[0]+1932/2197*h*K1[0]-7200/2197*h*K2[0]+7296/2197*h*K3[0])**2\
        *(secondVec[0]+1932/2197*h*M1[0]-7200/2197*h*M2[0]+7296/2197*h*M3[0])
    )
    #n=1,...,N-1
    for n in range(0,N):
        K4.append(
            -1j*g*np.abs(firstVec[n]+1932/2197*h*K1[n]-7200/2197*h*K2[n]+7296/2197*h*K3[n])**2\
            *(secondVec[n-1]+1932/2197*h*M1[n-1]-7200/2197*h*M2[n-1]+7296/2197*h*M3[n-1]
              +secondVec[n]+1932/2197*h*M1[n]-7200/2197*h*M2[n]+7296/2197*h*M3[n])
        )

    #n=0,...,N-2
    for n in  range(0,N-1):
        M4.append(
            -1j*g*np.abs(secondVec[n]+1932/2197*h*M1[n]-7200/2197*h*M2[n]+7296/2197*h*M3[n])**2\
            *(firstVec[n]+1932/2197*h*K1[n]-7200/2197*h*K2[n]+7296/2197*h*K3[n]
              +firstVec[n+1]+1932/2197*h*K1[n+1]-7200/2197*h*K2[n+1]+7296/2197*h*K3[n+1])
        )
    M4.append(
        -1j*g*np.abs(secondVec[N-1]+1932/2197*h*M1[N-1]-7200/2197*h*M2[N-1]+7296/2197*h*M3[N-1])**2\
        *(firstVec[N-1]+1932/2197*h*K1[N-1]-7200/2197*h*K2[N-1]+7296/2197*h*K3[N-1])
    )

    return K4, M4


def fifthVecRk5(firstVec,secondVec,g,h,K1,M1,K2,M2,K3,M3,K4,M4):
    '''

    :param firstVec: wvfcnt at even sites
    :param secondVec:  wvfcnt at odd sites
    :param g: nonlinearity
    :param h: time step
    :param K1: K1 coef vector
    :param M1: M1 coef vector
    :param K2: K2 coef vetcor
    :param M2: M2 coef vector
    :param K3: K3 coef vector
    :param M3: M3 coef vecror
    :param K4: K4 coef vector
    :param M4: M4 coef vector
    :return: K5, M5
    '''
    K5=[]
    M5=[]
    K5.append(
        -1j*g*np.abs(firstVec[0]+439/216*h*K1[0]-8*h*K2[0]+3680/513*h*K3[0]-845/4104*h*K4[0])**2\
        *(secondVec[0]+439/216*h*M1[0]-8*h*M2[0]+3680/513*h*M3[0]-854/4104*h*M4[0])
    )
    #n=1,...,N-1
    for n in range(1,N):
        K5.append(
            -1j*g*np.abs(firstVec[n]+439/216*h*K1[n]-8*h*K2[n]+3680/513*h*K3[n]-845/4104*h*K4[n])**2\
            *(secondVec[n-1]+439/216*h*M1[n-1]-8*h*M2[n-1]+3680/513*h*M3[n-1]-845/4104*h*M4[n-1]
              +secondVec[n]+439/216*h*M1[n]-8*h*M2[n]+3680/513*h*M3[n]-845/4104*h*M4[n])
        )
    #n=0,...,N-2
    for n in range(0,N-1):
        M5.append(
            -1j*g*np.abs(secondVec[n]+439/216*h*M1[n]-8*h*M2[n]+3680/513*h*M3[n]-845/4104*h*M4[n])**2\
            *(firstVec[n]+439/216*h*K1[n]-8*h*K2[n]+3680/513*h*K3[n]-845/4104*h*K4[n]
              +firstVec[n+1]+439/216*h*K1[n+1]-8*h*K2[n+1]+3680/513*h*K3[n+1]-845/4104*h*K4[n+1])
        )
    M5.append(
        -1j*g*np.abs(secondVec[N-1]+439/216*h*M1[N-1]-8*h*M2[N-1]+3680/513*h*M3[N-1]-845/4104*h*M4[N-1])**2\
        *(firstVec[N-1]+439/216*h*K1[N-1]-8*h*K2[N-1]+3680/513*h*K3[N-1]-845/4104*h*K4[N-1])
    )
    return K5, M5


