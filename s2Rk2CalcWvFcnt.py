from s2Rk2Fcnt import *
import inspect


gVals=[0]
outDir="./pump1/"
tStart = datetime.now()
for g in gVals:
    dataAll=[]
    dataAll.append(psi0)
    #init wvfcnt at even sites, odd sites
    a=[]
    b=[]
    for n in range(0,N):
        a.append(psi0[2*n])
        b.append(psi0[2*n+1])
    for q in range(0,Q):
        # time evolution starting from time step q=0,1,...,Q-1
        a,b=oneStepS2(a,b,q,g)
        psiCurr=[]
        for n in range(0,N):
            psiCurr.append(a[n])
            psiCurr.append(b[n])
        dataAll.append(psiCurr)

    xPos=[]
    for q in range(0,Q+1):
        psiTmp=dataAll[q]
        xMeanTmp=meanX(psiTmp)
        xPos.append(xMeanTmp)

    drift=[elem-xc for elem in xPos]
    tAll=[dt*q for q in range(0,Q+1)]
    plt.figure(figsize=(20,20))
    plt.plot(tAll,drift,color="black")
    plt.xlabel("time")
    plt.ylabel("ave position")
    plt.title("g = " + str(g) + ", width = " + str(sgm) + ", pumping = " + str(drift[-1] - drift[0]))
    plt.savefig(outDir + "g" + str(g) + "width" + str(sgm) + "position.png")
    plt.close()
    # write params info

    outTxt = outDir + "g" + str(g) + "width" + str(sgm) + "info.txt"

    fptr = open(outTxt, "w+")
    fptr.write("Total drift = " + str(drift[-1] - drift[0]) + "\n")
    fptr.write("g=" + str(g) + "\n")
    fptr.write("omega=" + str(omega) + "\n")
    fptr.write("omegaF=" + str(omegaF) + "\n")
    fptr.write("d0=" + str(d0) + "\n")
    fptr.write("D0=" + str(D0) + "\n")
    fptr.write("J=" + str(J) + "\n")
    fptr.write("k0=" + str(k0) + "\n")
    fptr.write(inspect.getsource(x))
    fptr.write(inspect.getsource(y))
    fptr.write(inspect.getsource(u))
    fptr.write(inspect.getsource(v))
    fptr.write(inspect.getsource(w))
    fptr.close()


tEnd = datetime.now()
print("computation time: ", tEnd - tStart)