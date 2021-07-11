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


def fifthVector(q,h,g,aVec,bVec,K1,M1,K2,M2,K3,M3,K4,M4):
    '''

    :param q: time step
    :param h: time duration
    :param g: nonlinearity
    :param aVec: input wvfcnt at even sites
    :param bVec: input wvfcnt at odd sites
    :param K1:  K1 coef vector
    :param M1: M1 coef vector
    :param K2: K2 coef vector
    :param M2: M2 coef vector
    :param K3: K3 coef vector
    :param M3: M3 coef vector
    :param K4: K4 coef vector
    :param M4: M4 coef vector
    :return: K5,M5
    '''
    K5=[]
    M5=[]
    timeQuad=q*h+(7-c)/14*h


    K5.append(
        -1j*v(timeQuad)*(bVec[0]+(9*c-21)/392*h*M1[0]+(8*c-56)/392*h*M2[0]
                         +(336-48*c)/392*h*M3[0]+(3*c-63)/392*h(M4[0]))\
        -1j*u(timeQuad)*(aVec[0]+(9*c-21)/392*h*K1[0]+(8*c-56)/392*h*K2[0]
                         +(336-48*c)/392*h*K3[0]+(3*c-63)/392*h*K4[0])\
        -1j*xmValsAll[0]*(aVec[0]+(9*c-21)/392*h*K1[0]+(8*c-56)/392*h*K2[0]
                         +(336-48*c)/392*h*K3[0]+(3*c-63)/392*h*K4[0])\
        -1j*g*np.abs(aVec[0]+(9*c-21)/392*h*K1[0]+(8*c-56)/392*h*K2[0]
                         +(336-48*c)/392*h*K3[0]+(3*c-63)/392*h*K4[0])**2
        *(bVec[0]+(9*c-21)/392*h*M1[0]+(8*c-56)/392*h*M2[0]
                         +(336-48*c)/392*h*M3[0]+(3*c-63)/392*h(M4[0]))
    )
    #n=1,...,N-1
    for n in range(1,N):
        K5.append(
            -1j*v(timeQuad)*(bVec[n]+(9*c-21)/392*h*M1[n]+(8*c-56)/392*h*M2[n]
                             +(336-48*c)/392*h*M3[n]+(3*c-63)/392*h*M4[n])\
            -1j*w(timeQuad)*(bVec[n-1]+(9*c-21)/392*h*M1[n-1]+(8*c-56)/392*h*M2[n-1]
                             +(336-48*c)/392*h*M3[n-1]+(3*c-63)/392*h*M4[n-1])\
            -1j*u(timeQuad)*(aVec[n]+(9*c-21)/392*h*K1[n]+(8*c-56)/392*h*K2[n]
                             +(336-48*c)/392*h*K3[n]+(3*c-63)/392*h*K4[n])\
            -1j*xmValsAll[n]*(aVec[n]+(9*c-21)/392*h*K1[n]+(8*c-56)/392*h*K2[n]
                             +(336-48*c)/392*h*K3[n]+(3*c-63)/392*h*K4[n])\
            -1j*g*np.abs(aVec[n]+(9*c-21)/392*h*K1[n]+(8*c-56)/392*h*K2[n]
                             +(336-48*c)/392*h*K3[n]+(3*c-63)/392*h*K4[n])**2\
            *(bVec[n-1]+(9*c-21)/392*h*M1[n-1]+(8*c-56)/392*h*M2[n-1]
                             +(336-48*c)/392*h*M3[n-1]+(3*c-63)/392*h*M4[n-1]
              +bVec[n]+(9*c-21)/392*h*M1[n]+(8*c-56)/392*h*M2[n]
                             +(336-48*c)/392*h*M3[n]+(3*c-63)/392*h*M4[n])
        )
    #n=0,...,N-2
    for n in range(0,N-1):
        M5.append(
            -1j*v(timeQuad)*(aVec[n]+(9*c-21)/392*h*K1[n]+(8*c-56)/392*h*K2[n]
                             +(336-48*c)/392*h*K3[n]+(3*c-63)/392*h*K4[n])\
            -1j*w(timeQuad)*(aVec[n+1]+(9*c-21)/392*h*K1[n+1]+(8*c-56)/392*h*K2[n+1]
                             +(336-48*c)/392*h*K3[n+1]+(3*c-63)/392*h*K4[n+1])\
            +1j*u(timeQuad)*(bVec[n]+(9*c-21)/392*h*M1[n]+(8*c-56)/392*h*M2[n]
                             +(336-48*c)/392*h*M3[n]+(3*c-63)/392*h*M4[n])\
            -1j*ymValsAll[n]*(bVec[n]+(9*c-21)/392*h*M1[n]+(8*c-56)/392*h*M2[n]
                             +(336-48*c)/392*h*M3[n]+(3*c-63)/392*h*M4[n])\
            -1j*g*np.abs(bVec[n]+(9*c-21)/392*h*M1[n]+(8*c-56)/392*h*M2[n]
                             +(336-48*c)/392*h*M3[n]+(3*c-63)/392*h*M4[n])**2\
            *(aVec[n]+(9*c-21)/392*h*K1[n]+(8*c-56)/392*h*K2[n]
                             +(336-48*c)/392*h*K3[n]+(3*c-63)/392*h*K4[n]
              +aVec[n+1]+(9*c-21)/392*h*K1[n+1]+(8*c-56)/392*h*K2[n+1]
                             +(336-48*c)/392*h*K3[n+1]+(3*c-63)/392*h*K4[n+1])
        )

    M5.append(
        -1j*v(timeQuad)*(aVec[N-1]+(9*c-21)/392*h*K1[N-1]+(8*c-56)/392*h*K2[N-1]
                         +(336-48*c)/392*h*K3[N-1]+(3*c-63)/392*h*K4[N-1])\
        +1j*u(timeQuad)*(bVec[N-1]+(9*c-21)/392*h*M1[N-1]+(8*c-56)/392*h*M2[N-1]
                         +(336-48*c)/392*h*M3[N-1]+(3*c-63)/392*h*M4[N-1])\
        -1j*ymValsAll[N-1]*(bVec[N-1]+(9*c-21)/392*h*M1[N-1]+(8*c-56)/392*h*M2[N-1]
                         +(336-48*c)/392*h*M3[N-1]+(3*c-63)/392*h*M4[N-1])\
        -1j*g*np.abs(bVec[N-1]+(9*c-21)/392*h*M1[N-1]+(8*c-56)/392*h*M2[N-1]
                         +(336-48*c)/392*h*M3[N-1]+(3*c-63)/392*h*M4[N-1])**2\
        *(aVec[N-1]+(9*c-21)/392*h*K1[N-1]+(8*c-56)/392*h*K2[N-1]
                         +(336-48*c)/392*h*K3[N-1]+(3*c-63)/392*h*K4[N-1])
    )

    return K5, M5

def sixthVector(q,h,g,aVec,bVec,K1,M1,K2,M2,K3,M3,K4,M4,K5,M5):
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
    :param K4: K4 coef vector
    :param M4: M4 coef vector
    :param K5: K5 coef vector
    :param M5: M5 coef vector
    :return: K6, M6
    '''
    K6=[]
    M6=[]
    timeQuad=q*h+(7+c)/14*h

    K6.append(
        -1j*v(timeQuad)*(bVec[0]-(1155+255*c)/1960*h*M1[0]-(280+40*c)/1960*h*M2[0]
                         -320*c/1960*h*M3[0]+(63+363*c)/1960*h*M4[0]+(2352+392*c)/1960*h*M5[0])\
        -1j*u(timeQuad)*(aVec[0]-(1155+255*c)/1960*h*K1[0]-(280+40*c)/1960*h*K2[0]
                         -320*c/1960*h*K3[0]+(63+363*c)/1960*h*K4[0]+(2352+392*c)/1960*h*K5[0])\
        -1j*xmValsAll[0]*(aVec[0]-(1155+255*c)/1960*h*K1[0]-(280+40*c)/1960*h*K2[0]
                         -320*c/1960*h*K3[0]+(63+363*c)/1960*h*K4[0]+(2352+392*c)/1960*h*K5[0])\
        -1j*g*np.abs(aVec[0]-(1155+255*c)/1960*h*K1[0]-(280+40*c)/1960*h*K2[0]
                         -320*c/1960*h*K3[0]+(63+363*c)/1960*h*K4[0]+(2352+392*c)/1960*h*K5[0])**2\
        *(bVec[0]-(1155+255*c)/1960*h*M1[0]-(280+40*c)/1960*h*M2[0]
                         -320*c/1960*h*M3[0]+(63+363*c)/1960*h*M4[0]+(2352+392*c)/1960*h*M5[0])
    )
    #n=1,...,N-1
    for n in range(1,N):
        K6.append(
            -1j*v(timeQuad)*(bVec[n]-(1155+255*c)/1960*h*M1[n]-(280+40*c)/1960*h*M2[n]
                             -320*c/1960*h*M3[n]+(63+363*c)/1960*h*M4[n]+(2352+392*c)/1960*h*M5[n])\
            -1j*w(timeQuad)*(bVec[n-1]-(1155+255*c)/1960*h*M1[n-1]-(280+40*c)/1960*h*M2[n-1]
                             -320*c/1960*h*M3[n-1]+(63+363*c)/1960*h*M4[n-1]+(2352+392*c)/1960*h*M5[n-1])\
            -1j*u(timeQuad)*(aVec[n]-(1155+255*c)/1960*h*K1[n]-(280+40*c)/1960*h*K2[n]
                             -320*c/1960*K3[n]+(63+363*c)/1960*h*K4[n]+(2352+392*c)/1960*h*K5[n])\
            -1j*xmValsAll[n]*(aVec[n]-(1155+255*c)/1960*h*K1[n]-(280+40*c)/1960*h*K2[n]
                             -320*c/1960*K3[n]+(63+363*c)/1960*h*K4[n]+(2352+392*c)/1960*h*K5[n])\
            -1j*g*np.abs(aVec[n]-(1155+255*c)/1960*h*K1[n]-(280+40*c)/1960*h*K2[n]
                             -320*c/1960*K3[n]+(63+363*c)/1960*h*K4[n]+(2352+392*c)/1960*h*K5[n])**2\
            *(bVec[n-1]-(1155+255*c)/1960*h*M1[n-1]-(280+40*c)/1960*h*M2[n-1]
                             -320*c/1960*h*M3[n-1]+(63+363*c)/1960*h*M4[n-1]+(2352+392*c)/1960*h*M5[n-1]
              +bVec[n]-(1155+255*c)/1960*h*M1[n]-(280+40*c)/1960*h*M2[n]
                             -320*c/1960*h*M3[n]+(63+363*c)/1960*h*M4[n]+(2352+392*c)/1960*h*M5[n])
        )
    #n=0,...,N-2
    for n in range(0,N-1):
        M6.append(
            -1j*v(timeQuad)*(aVec[n]-(1155+255*c)/1960*h*K1[n]-(280+40*c)/1960*h*K2[n]
                             -320*c/1960*h*K3[n]+(63+363*c)/1960*h*K4[n]+(2352+392*c)/1960*h*K5[n])\
            -1j*w(timeQuad)*(aVec[n+1]-(1155+255*c)/1960*h*K1[n+1]-(280+40*c)/1960*h*K2[n+1]
                             -320*c/1960*h*K3[n+1]+(63+363*c)/1960*h*K4[n+1]+(2352+392*c)/1960*h*K5[n+1])\
            +1j*u(timeQuad)*(bVec[n]-(1155+255*c)/1960*h*M1[n]-(280+40*c)/1960*h*M2[n]
                             -320*c/1960*h*M3[n]+(63+363*c)/1960*h*M4[n]+(2352+392*c)/1960*h*M5[n])\
            -1j*ymValsAll[n]*(bVec[n]-(1155+255*c)/1960*h*M1[n]-(280+40*c)/1960*h*M2[n]
                             -320*c/1960*h*M3[n]+(63+363*c)/1960*h*M4[n]+(2352+392*c)/1960*h*M5[n])\
            -1j*g*np.abs(bVec[n]-(1155+255*c)/1960*h*M1[n]-(280+40*c)/1960*h*M2[n]
                             -320*c/1960*h*M3[n]+(63+363*c)/1960*h*M4[n]+(2352+392*c)/1960*h*M5[n])**2\
            *(aVec[n]-(1155+255*c)/1960*h*K1[n]-(280+40*c)/1960*h*K2[n]
                             -320*c/1960*h*K3[n]+(63+363*c)/1960*h*K4[n]+(2352+392*c)/1960*h*K5[n]
              +aVec[n+1]-(1155+255*c)/1960*h*K1[n+1]-(280+40*c)/1960*h*K2[n+1]
                             -320*c/1960*h*K3[n+1]+(63+363*c)/1960*h*K4[n+1]+(2352+392*c)/1960*h*K5[n+1])
        )

    M6.append(
        -1j*v(timeQuad)*(aVec[N-1]-(1155+255*c)/1960*h*K1[N-1]-(280+40*c)/1960*h*K2[N-1]
                         -320*c/1960*h*K3[N-1]+(63+363*c)/1960*h*K4[N-1]+(2352+392*c)/1960*h*K5[N-1])\
        +1j*u(timeQuad)*(bVec[N-1]-(1155+255*c)/1960*h*M1[N-1]-(280+40*c)/1960*h*M2[N-1]
                         -320*c/1960*h*M3[N-1]+(63+363*c)/1960*h*M4[N-1]+(2352+392*c)/1960*h*M5[N-1])\
        -1j*ymValsAll[N-1]*(bVec[N-1]-(1155+255*c)/1960*h*M1[N-1]-(280+40*c)/1960*h*M2[N-1]
                         -320*c/1960*h*M3[N-1]+(63+363*c)/1960*h*M4[N-1]+(2352+392*c)/1960*h*M5[N-1])\
        -1j*g*np.abs(bVec[N-1]-(1155+255*c)/1960*h*M1[N-1]-(280+40*c)/1960*h*M2[N-1]
                         -320*c/1960*h*M3[N-1]+(63+363*c)/1960*h*M4[N-1]+(2352+392*c)/1960*h*M5[N-1])**2\
        *(aVec[N-1]-(1155+255*c)/1960*h*K1[N-1]-(280+40*c)/1960*h*K2[N-1]
                         -320*c/1960*h*K3[N-1]+(63+363*c)/1960*h*K4[N-1]+(2352+392*c)/1960*h*K5[N-1])
    )
    return K6,M6


def seventhVector(q,h,g,aVec,bVec,K1,M1,K2,M2,K3,M3,K4,M4,K5,M5,K6,M6):
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
    :param K4: K4 coef vector
    :param M4: M4 coef vector
    :param K5: K5 coef vector
    :param M5: M5 coef vector
    :param K6: K6 coef vector
    :param M6: M6 coef vector
    :return: K6, M6
    '''
    K7=[]
    M7=[]
    timeQuad=q*h+h

    K7.append(
        -1j*v(timeQuad)*(bVec[0]+(330+105*c)/180*h*M1[0]+120*h*M2[0]+(280*c-200)/180*h*M3[0]
                         +(126-189*c)/180*h*M4[0]-(686+126*c)/180*h*M5[0]+(490-70*c)/180*h*M6[0])\
        -1j*u(timeQuad)*(aVec[0]+(330+105*c)/180*h*K1[0]+120*h*K2[0]+(280*c-200)/180*h*K3[0]
                         +(126-189*c)/180*h*K4[0]-(686+126*c)/180*h*K5[0]+(490-70*c)/180*h*K6[0])\
        -1j*xmValsAll[0]*(aVec[0]+(330+105*c)/180*h*K1[0]+120*h*K2[0]+(280*c-200)/180*h*K3[0]
                         +(126-189*c)/180*h*K4[0]-(686+126*c)/180*h*K5[0]+(490-70*c)/180*h*K6[0])\
        -1j*g*np.abs(aVec[0]+(330+105*c)/180*h*K1[0]+120*h*K2[0]+(280*c-200)/180*h*K3[0]
                         +(126-189*c)/180*h*K4[0]-(686+126*c)/180*h*K5[0]+(490-70*c)/180*h*K6[0])**2\
        *(bVec[0]+(330+105*c)/180*h*M1[0]+120*h*M2[0]+(280*c-200)/180*h*M3[0]
                         +(126-189*c)/180*h*M4[0]-(686+126*c)/180*h*M5[0]+(490-70*c)/180*h*M6[0])
    )

    #n=1,...,N-1
    for n in range(1,N):
        K7.append(
            -1j*v(timeQuad)*(bVec[n]+(330+105*c)/180*h*M1[n]+120*h*M2[n]+(280*c-200)/180*h*M3[n]
                             +(126-189*c)/180*h*M4[n]-(686+126*c)/180*h*M5[n]+(490-70*c)/180*h*M6[n])\
            -1j*w(timeQuad)*(bVec[n-1]+(330+105*c)/180*h*M1[n-1]+120*h*M2[n-1]+(280*c-200)/180*h*M3[n-1]
                             +(126-189*c)/180*h*M4[n-1]-(686+126*c)/180*h*M5[n-1]+(490-70*c)/180*h*M6[n-1])\
            -1j*u(timeQuad)*(aVec[n]+(330+105*c)/180*h*K1[n]+120*h*K2[n]+(280*c-200)/180*h*K3[n]
                             +(126-189*c)/180*h*K4[n]-(686+126*c)/180*h*K5[n]+(490-70*c)/180*h*K6[n])\
            -1j*xmValsAll[n]*(aVec[n]+(330+105*c)/180*h*K1[n]+120*h*K2[n]+(280*c-200)/180*h*K3[n]
                             +(126-189*c)/180*h*K4[n]-(686+126*c)/180*h*K5[n]+(490-70*c)/180*h*K6[n])\
            -1j*g*np.abs(aVec[n]+(330+105*c)/180*h*K1[n]+120*h*K2[n]+(280*c-200)/180*h*K3[n]
                             +(126-189*c)/180*h*K4[n]-(686+126*c)/180*h*K5[n]+(490-70*c)/180*h*K6[n])**2\
            *(bVec[n-1]+(330+105*c)/180*h*M1[n-1]+120*h*M2[n-1]+(280*c-200)/180*h*M3[n-1]
                             +(126-189*c)/180*h*M4[n-1]-(686+126*c)/180*h*M5[n-1]+(490-70*c)/180*h*M6[n-1]
              +bVec[n]+(330+105*c)/180*h*M1[n]+120*h*M2[n]+(280*c-200)/180*h*M3[n]
                             +(126-189*c)/180*h*M4[n]-(686+126*c)/180*h*M5[n]+(490-70*c)/180*h*M6[n])
        )
    #n=0,...,N-2
    for n in range(0,N-1):
        M7.append(
            -1j*v(timeQuad)*(aVec[n]+(330+105*c)/180*h*K1[n]+120*h*K2[n]+(280*c-200)/180*h*K3[n]
                             +(126-189*c)/180*h*K4[n]-(686+126*c)/180*h*K5[n]+(490-70*c)/180*h*K6[n])\
            -1j*w(timeQuad)*(aVec[n+1]+(330+105*c)/180*h*K1[n+1]+120*h*K2[n+1]+(280*c-200)/180*h*K3[n+1]
                             +(126-189*c)/180*h*K4[n+1]-(686+126*c)/180*h*K5[n+1]+(490-70*c)/180*h*K6[n+1])\
            +1j*u(timeQuad)*(bVec[n]+(330+105*c)/180*h*M1[n]+120*h*M2[n]+(280*c-200)/180*h*M3[n]
                              +(126-189*c)/180*h*M4[n]-(686+126*c)/180*h*M5[n]+(490-70*c)/180*h*M6[n])\
            -1j*ymValsAll[n]*(bVec[n]+(330+105*c)/180*h*M1[n]+120*h*M2[n]+(280*c-200)/180*h*M3[n]
                              +(126-189*c)/180*h*M4[n]-(686+126*c)/180*h*M5[n]+(490-70*c)/180*h*M6[n])\
            -1j*g*np.abs(bVec[n]+(330+105*c)/180*h*M1[n]+120*h*M2[n]+(280*c-200)/180*h*M3[n]
                              +(126-189*c)/180*h*M4[n]-(686+126*c)/180*h*M5[n]+(490-70*c)/180*h*M6[n])**2\
            *(aVec[n]+(330+105*c)/180*h*K1[n]+120*h*K2[n]+(280*c-200)/180*h*K3[n]
                             +(126-189*c)/180*h*K4[n]-(686+126*c)/180*h*K5[n]+(490-70*c)/180*h*K6[n]
              +aVec[n+1]+(330+105*c)/180*h*K1[n+1]+120*h*K2[n+1]+(280*c-200)/180*h*K3[n+1]
                             +(126-189*c)/180*h*K4[n+1]-(686+126*c)/180*h*K5[n+1]+(490-70*c)/180*h*K6[n+1])
        )
    M7.append(
        -1j*v(timeQuad)*(aVec[N-1]+(330+105*c)/180*h*K1[N-1]+120*h*K2[N-1]+(280*c-200)/180*h*K3[N-1]
                         +(126-189*c)/180*h*K4[N-1]-(686+126*c)/180*h*K5[N-1]+(490-70*c)/180*h*K6[N-1])\
        +1j*u(timeQuad)*(bVec[N-1]+(330+105*c)/180*h*M1[N-1]+120*h*M2[N-1]+(280*c-200)/180*h*M3[N-1]
                              +(126-189*c)/180*h*M4[N-1]-(686+126*c)/180*h*M5[N-1]+(490-70*c)/180*h*M6[N-1])\
        -1j*ymValsAll[N-1]*(bVec[N-1]+(330+105*c)/180*h*M1[N-1]+120*h*M2[N-1]+(280*c-200)/180*h*M3[N-1]
                              +(126-189*c)/180*h*M4[N-1]-(686+126*c)/180*h*M5[N-1]+(490-70*c)/180*h*M6[N-1])\
        -1j*g*np.abs(bVec[N-1]+(330+105*c)/180*h*M1[N-1]+120*h*M2[N-1]+(280*c-200)/180*h*M3[N-1]
                              +(126-189*c)/180*h*M4[N-1]-(686+126*c)/180*h*M5[N-1]+(490-70*c)/180*h*M6[N-1])**2\
        *(aVec[N-1]+(330+105*c)/180*h*K1[N-1]+120*h*K2[N-1]+(280*c-200)/180*h*K3[N-1]
                         +(126-189*c)/180*h*K4[N-1]-(686+126*c)/180*h*K5[N-1]+(490-70*c)/180*h*K6[N-1])


    )
    return  K7, M7



def oneStepRk6(q,g,aVec,bVec):
    K1, M1=firstVector(q,dt,g,aVec,bVec)
    K2,M2=secondVector(q,dt,g,aVec,bVec,K1,M1)
    K3,M3=thirdVector(q,dt,g,aVec,bVec,K1,M1,K2,M2)
    K4,M4=fourthVector(q,dt,g,aVec,bVec,K1,M1,K2,M2,K3,M3)
    K5,M5=fifthVector(q,dt,g,aVec,bVec,K1,M1,K2,M2,K3,M3,K4,M4)
    K6,M6=sixthVector(q,dt,g,aVec,bVec,K1,M1,K2,M2,K3,M3,K4,M4,K5,M5)
    K7,M7=seventhVector(q,dt,g,aVec,bVec,K1,M1,K2,M2,K3,M3,K4,M4,K5,M5,K6,M6)
    aNext=[]
    bNext=[]
    for n in range(0,N):
        aNext.append(
            aVec[n]+1/180*dt*(9*K1[n]+64*K3[n]+49*K5[n]+49*K6[n]+9*K7[n])
        )
    for n in range(0,N):
        bNext.append(
            bVec[n]+1/180*dt*(9*M1[n]+64*M3[n]+49*M5[n]+49*M6[n]+9*M7[n])
        )

    return aNext,bNext