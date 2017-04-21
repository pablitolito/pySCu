#All of this cade coming from PmagPy software libraries (Tauxe et al. 2016, G3, DOI: 10.1002/2016GC006307)

import  numpy,string,sys
from numpy import random
import numpy.linalg
import exceptions
import os
#import check_updates
#import scipy
#from scipy import array,sqrt,mean

#check_updates.main() # check for updates

def flip(D):
    """
     flip reverse mode
    """
    ppars=doprinc(D) # get principle direction
    D1,D2=[],[]
    for rec in D:
        ang=angle([rec[0],rec[1]],[ppars['dec'],ppars['inc']])
        if ang>90.:
            d,i=(rec[0]-180.)%360.,-rec[1]
            D2.append([d,i,1.])
        else:
            D1.append([rec[0],rec[1],1.])
    return D1,D2

def dotilt(dec,inc,bed_az,bed_dip):
    """
    does a tilt correction on dec,inc using bedding dip direction bed_az and dip bed_dip.  called with syntax:  dotilt(dec,inc,bed_az,bed_dip).
    """
    rad=numpy.pi/180. # converts from degrees to radians
    X=dir2cart([dec,inc,1.]) # get cartesian coordinates of dec,inc
# get some sines and cosines of new coordinate system
    sa,ca= -numpy.sin(bed_az*rad),numpy.cos(bed_az*rad) 
    cdp,sdp= numpy.cos(bed_dip*rad),numpy.sin(bed_dip*rad) 
# do the rotation
    xc=X[0]*(sa*sa+ca*ca*cdp)+X[1]*(ca*sa*(1.-cdp))+X[2]*sdp*ca
    yc=X[0]*ca*sa*(1.-cdp)+X[1]*(ca*ca+sa*sa*cdp)-X[2]*sa*sdp
    zc=X[0]*ca*sdp-X[1]*sdp*sa-X[2]*cdp
# convert back to direction:
    Dir=cart2dir([xc,yc,-zc])
    return Dir[0],Dir[1] # return declination, inclination of rotated direction

def dotilt_V(input):
    """
    does a tilt correction on dec,inc using bedding dip direction bed_az and dip bed_dip
    """
    input=input.transpose() 
    dec, inc, bed_az, bed_dip =input[0],input[1],input[2],input[3]  # unpack input array into separate arrays
    rad=numpy.pi/180. # convert to radians
    Dir=numpy.array([dec,inc]).transpose()
    X=dir2cart(Dir).transpose() # get cartesian coordinates
    N=numpy.size(dec)

# get some sines and cosines of new coordinate system
    sa,ca= -numpy.sin(bed_az*rad),numpy.cos(bed_az*rad) 
    cdp,sdp= numpy.cos(bed_dip*rad),numpy.sin(bed_dip*rad) 
# do the rotation
    xc=X[0]*(sa*sa+ca*ca*cdp)+X[1]*(ca*sa*(1.-cdp))+X[2]*sdp*ca
    yc=X[0]*ca*sa*(1.-cdp)+X[1]*(ca*ca+sa*sa*cdp)-X[2]*sa*sdp
    zc=X[0]*ca*sdp-X[1]*sdp*sa-X[2]*cdp
# convert back to direction:
    cart=numpy.array([xc,yc,-zc]).transpose()
    Dir=cart2dir(cart).transpose()
    return Dir[0],Dir[1] # return declination, inclination arrays of rotated direction

def dogeo(dec,inc,az,pl):
    """
    rotates dec,in into geographic coordinates using az,pl as azimuth and plunge of X direction
    """
    A1,A2,A3=[],[],[] # set up lists for rotation vector
    Dir=[dec,inc,1.] # put dec inc in direction list and set  length to unity
    X=dir2cart(Dir) # get cartesian coordinates
#
#   set up rotation matrix
#
    A1=dir2cart([az,pl,1.])
    A2=dir2cart([az+90.,0,1.])
    A3=dir2cart([az-180.,90.-pl,1.])
#
# do rotation
#
    xp=A1[0]*X[0]+A2[0]*X[1]+A3[0]*X[2]
    yp=A1[1]*X[0]+A2[1]*X[1]+A3[1]*X[2]
    zp=A1[2]*X[0]+A2[2]*X[1]+A3[2]*X[2]
#
# transform back to dec,inc
#
    Dir_geo=cart2dir([xp,yp,zp])
    return Dir_geo[0],Dir_geo[1]    # send back declination and inclination

def dogeo_V(input):
    """
    rotates dec,in into geographic coordinates using az,pl as azimuth and plunge of X direction
    handles  array for  input 
    """
    input=input.transpose() 
    dec, inc, az, pl =input[0],input[1],input[2],input[3]  # unpack input array into separate arrays
    Dir=numpy.array([dec,inc]).transpose()
    X=dir2cart(Dir).transpose() # get cartesian coordinates
    N=numpy.size(dec)
    A1=dir2cart(numpy.array([az,pl,numpy.ones(N)]).transpose()).transpose()
    A2=dir2cart(numpy.array([az+90.,numpy.zeros(N),numpy.ones(N)]).transpose()).transpose()
    A3=dir2cart(numpy.array([az-180.,90.-pl,numpy.ones(N)]).transpose()).transpose()

# do rotation
#
    xp=A1[0]*X[0]+A2[0]*X[1]+A3[0]*X[2]
    yp=A1[1]*X[0]+A2[1]*X[1]+A3[1]*X[2]
    zp=A1[2]*X[0]+A2[2]*X[1]+A3[2]*X[2]
    cart=numpy.array([xp,yp,zp]).transpose()
#
# transform back to dec,inc
#
    Dir_geo=cart2dir(cart).transpose()
    return Dir_geo[0],Dir_geo[1]    # send back declination and inclination arrays

def dodirot(D,I,Dbar,Ibar):
    """
    This function is called by dodirot(D,I,Dbar,Ibar) where D=declination, I = inclination and Dbar/Ibar are the desired mean direction.  It returns the rotated Dec/Inc pair.
    """
    d,irot=dogeo(D,I,Dbar,90.-Ibar)
    drot=d-180.
#    drot,irot=dogeo(D,I,Dbar,Ibar)
    if drot<360.:drot=drot+360.
    if drot>360.:drot=drot-360.
    return drot,irot

def angle(D1,D2):
    """
    call to angle(D1,D2) returns array of angles between lists of two directions D1,D2 where D1 is for example, [[Dec1,Inc1],[Dec2,Inc2],etc.]
    """
    D1=numpy.array(D1)
    if len(D1.shape)>1:
        D1=D1[:,0:2] # strip off intensity
    else: D1=D1[:2]
    D2=numpy.array(D2)
    if len(D2.shape)>1:
        D2=D2[:,0:2] # strip off intensity
    else: D2=D2[:2]
    X1=dir2cart(D1) # convert to cartesian from polar
    X2=dir2cart(D2)
    angles=[] # set up a list for angles
    for k in range(X1.shape[0]): # single vector
        angle= numpy.arccos(numpy.dot(X1[k],X2[k]))*180./numpy.pi # take the dot product
        angle=angle%360.
        angles.append(angle)
    return numpy.array(angles)

def cart2dir(cart):
    """
    converts a direction to cartesian coordinates
    """
    cart=numpy.array(cart)
    rad=numpy.pi/180. # constant to convert degrees to radians
    if len(cart.shape)>1:
        Xs,Ys,Zs=cart[:,0],cart[:,1],cart[:,2]
    else: #single vector
        Xs,Ys,Zs=cart[0],cart[1],cart[2]
    Rs=numpy.sqrt(Xs**2+Ys**2+Zs**2) # calculate resultant vector length
    Decs=(numpy.arctan2(Ys,Xs)/rad)%360. # calculate declination taking care of correct quadrants (arctan2) and making modulo 360.
    try:
        Incs=numpy.arcsin(Zs/Rs)/rad # calculate inclination (converting to degrees) # 
    except:
        print 'trouble in cart2dir' # most likely division by zero somewhere
        return numpy.zeros(3)
        
    return numpy.array([Decs,Incs,Rs]).transpose() # return the directions list

def tauV(T):
    """
    gets the eigenvalues (tau) and eigenvectors (V) from matrix T
    """
    t,V,tr=[],[],0.
    ind1,ind2,ind3=0,1,2
    evalues,evectmps=numpy.linalg.eig(T)
    evectors=numpy.transpose(evectmps)  # to make compatible with Numeric convention
    for tau in evalues:
        tr+=tau
    if tr!=0:
        for i in range(3):
            evalues[i]=evalues[i]/tr
    else:
        return t,V
# sort evalues,evectors
    t1,t2,t3=0.,0.,1.
    for k in range(3):
        if evalues[k] > t1: 
            t1,ind1=evalues[k],k 
        if evalues[k] < t3: 
            t3,ind3=evalues[k],k 
    for k in range(3):
        if evalues[k] != t1 and evalues[k] != t3: 
            t2,ind2=evalues[k],k
    V.append(evectors[ind1])
    V.append(evectors[ind2])
    V.append(evectors[ind3])
    t.append(t1)
    t.append(t2)
    t.append(t3)
    return t,V

def Tmatrix(X):
    """
    gets the orientation matrix (T) from data in X
    """
    T=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
    for row in X:
        for k in range(3):
            for l in range(3):
                T[k][l] += row[k]*row[l]
    return T

def dir2cart(d):
   # converts list or array of vector directions, in degrees, to array of cartesian coordinates, in x,y,z
    ints=numpy.ones(len(d)).transpose() # get an array of ones to plug into dec,inc pairs
    d=numpy.array(d)
    rad=numpy.pi/180.
    if len(d.shape)>1: # array of vectors
        decs,incs=d[:,0]*rad,d[:,1]*rad
        if d.shape[1]==3: ints=d[:,2] # take the given lengths
    else: # single vector
        decs,incs=numpy.array(d[0])*rad,numpy.array(d[1])*rad
        if len(d)==3: 
            ints=numpy.array(d[2])
        else:
            ints=numpy.array([1.])
    cart= numpy.array([ints*numpy.cos(decs)*numpy.cos(incs),ints*numpy.sin(decs)*numpy.cos(incs),ints*numpy.sin(incs)]).transpose()
    return cart

def circ(dec,dip,alpha):
    """
    function to calculate points on an circle about dec,dip with angle alpha
    """
    rad=numpy.pi/180.
    D_out,I_out=[],[]
    dec,dip,alpha=dec*rad ,dip*rad,alpha*rad
    dec1=dec+numpy.pi/2.
    isign=1
    if dip!=0: isign=(abs(dip)/dip)
    dip1=(dip-isign*(numpy.pi/2.))
    t=[[0,0,0],[0,0,0],[0,0,0]]
    v=[0,0,0]
    t[0][2]=numpy.cos(dec)*numpy.cos(dip)
    t[1][2]=numpy.sin(dec)*numpy.cos(dip)
    t[2][2]=numpy.sin(dip)
    t[0][1]=numpy.cos(dec)*numpy.cos(dip1)
    t[1][1]=numpy.sin(dec)*numpy.cos(dip1)
    t[2][1]=numpy.sin(dip1)
    t[0][0]=numpy.cos(dec1)
    t[1][0]=numpy.sin(dec1)
    t[2][0]=0   
    for i in range(101): 
        psi=float(i)*numpy.pi/50. 
        v[0]=numpy.sin(alpha)*numpy.cos(psi) 
        v[1]=numpy.sin(alpha)*numpy.sin(psi) 
        v[2]=numpy.sqrt(abs(1.-v[0]**2 - v[1]**2))
        elli=[0,0,0]
        for j in range(3):
            for k in range(3):
                elli[j]=elli[j] + t[j][k]*v[k] 
        Dir=cart2dir(elli)
        D_out.append(Dir[0])
        I_out.append(Dir[1])
    return D_out,I_out

def fisher_mean(data):
    """
    call to fisher_mean(data) calculates fisher statistics for data, which is a list of [dec,inc] pairs.  
    """

    R,Xbar,X,fpars=0,[0,0,0],[],{}
    N=len(data)
    if N <2:
       return fpars
    X=dir2cart(data)
    for i in range(len(X)):
        for c in range(3):
           Xbar[c]+=X[i][c]
    for c in range(3):
        R+=Xbar[c]**2
    R=numpy.sqrt(R)
    for c in range(3):
        Xbar[c]=Xbar[c]/R    
    dir=cart2dir(Xbar)
    fpars["dec"]=dir[0]
    fpars["inc"]=dir[1]
    fpars["n"]=N
    fpars["r"]=R
    if N!=R:
        k=(N-1.)/(N-R)
        fpars["k"]=k
        csd=81./numpy.sqrt(k)
    else:
        fpars['k']='inf'
        csd=0.
    b=20.**(1./(N-1.)) -1
    a=1-b*(N-R)/R
    if a<-1:a=-1
    a95=numpy.arccos(a)*180./numpy.pi
    fpars["alpha95"]=a95
    fpars["csd"]=csd
    if a<0: fpars["alpha95"] = 180.0
    return fpars
 
def gausspars(data):
    """
    calculates gaussian statistics for data
    """
    N,mean,d=len(data),0.,0.
    if N<1: return "",""
    if N==1: return data[0],0
    for j in range(N):
       mean+=data[j]/float(N)
    for j in range(N):
       d+=(data[j]-mean)**2 
    stdev=numpy.sqrt(d*(1./(float(N-1))))
    return mean,stdev

def weighted_mean(data):
    """
    calculates weighted mean of data
    """
    W,N,mean,d=0,len(data),0,0
    if N<1: return "",""
    if N==1: return data[0][0],0
    for x in data:
       W+=x[1] # sum of the weights
    for x in data:
       mean+=(float(x[1])*float(x[0]))/float(W)
    for x in data:
       d+=(float(x[1])/float(W))*(float(x[0])-mean)**2 
    stdev=numpy.sqrt(d*(1./(float(N-1))))
    return mean,stdev

def fisher_by_pol(data):
    """
    input:    as in dolnp (list of dictionaries with 'dec' and 'inc')
    description: do fisher mean after splitting data into two polarity domains.
    output: three dictionaries:
        'A'= polarity 'A'
        'B = polarity 'B'
        'ALL'= switching polarity of 'B' directions, and calculate fisher mean of all data     
    code modified from eqarea_ell.py b rshaar 1/23/2014
    """
    FisherByPoles={}
    DIblock,nameblock,locblock=[],[],[]
    for rec in data:
        if 'dec' in rec.keys() and 'inc' in rec.keys():
            DIblock.append([float(rec["dec"]),float(rec["inc"])]) # collect data for fisher calculation
        else:
            continue
        if 'name' in rec.keys():
            nameblock.append(rec['name'])
        else:
            nameblock.append("")    
        if 'loc' in rec.keys():
            locblock.append(rec['loc'])
        else:
            locblock.append("")
            
    ppars=doprinc(array(DIblock)) # get principal directions  
    reference_DI=[ppars['dec'],ppars['inc']] # choose the northerly declination principe component ("normal") 
    if reference_DI[0]>90 and reference_DI[0]<270: # make reference direction in northern hemisphere
        reference_DI[0]=(reference_DI[0]+180.)%360
        reference_DI[1]=reference_DI[1]*-1.
    nDIs,rDIs,all_DI,npars,rpars=[],[],[],[],[]
    nlist,rlist,alllist="","",""
    nloclist,rloclist,allloclist="","",""
    for k in range(len(DIblock)):            
        if angle([DIblock[k][0],DIblock[k][1]],reference_DI) > 90.:
            rDIs.append(DIblock[k])
            rlist=rlist+":"+nameblock[k]
            if locblock[k] not in rloclist:rloclist=rloclist+":"+locblock[k]
            all_DI.append( [(DIblock[k][0]+180.)%360.,-1.*DIblock[k][1]])
            alllist=alllist+":"+nameblock[k]
            if locblock[k] not in allloclist:allloclist=allloclist+":"+locblock[k]
        else:
            nDIs.append(DIblock[k])
            nlist=nlist+":"+nameblock[k]
            if locblock[k] not in nloclist:nloclist=nloclist+":"+locblock[k]
            all_DI.append(DIblock[k])
            alllist=alllist+":"+nameblock[k]
            if locblock[k] not in allloclist:allloclist=allloclist+":"+locblock[k]
            
    for mode in ['A','B','All']:
        if mode=='A' and len(nDIs)>2:
            fpars=fisher_mean(nDIs)
            fpars['sites']=nlist.strip(':')
            fpars['locs']=nloclist.strip(':')
            FisherByPoles[mode]=fpars
        elif mode=='B' and len(rDIs)>2:              
            fpars=fisher_mean(rDIs)
            fpars['sites']=rlist.strip(':')
            fpars['locs']=rloclist.strip(':')
            FisherByPoles[mode]=fpars
        elif mode=='All' and len(all_DI)>2:           
            fpars=fisher_mean(all_DI)
            fpars['sites']=alllist.strip(':')
            fpars['locs']=allloclist.strip(':')
            FisherByPoles[mode]=fpars
    return FisherByPoles       
    
def dolnp(data,direction_type_key):
    """
    returns fisher mean, a95 for data  using method of mcfadden and mcelhinny '88 for lines and planes
    """
    if "tilt_correction" in data[0].keys(): 
        tc=data[0]["tilt_correction"]
    else:
        tc='-1'
    n_lines,n_planes=0,0
    X,L,fdata,dirV=[],[],[],[0,0,0]
    E=[0,0,0]
    fpars={}
#
# sort data  into lines and planes and collect cartesian coordinates
    for rec in data:
        cart=dir2cart([rec["dec"],rec["inc"]])[0]
        if direction_type_key in rec.keys() and rec[direction_type_key]=='p': # this is a pole to a plane
            n_planes+=1
            L.append(cart) # this is the "EL, EM, EN" array of MM88
        else: # this is a line
            n_lines+=1
            fdata.append([rec["dec"],rec["inc"],1.]) # collect data for fisher calculation
            X.append(cart)
            E[0]+=cart[0] 
            E[1]+=cart[1] 
            E[2]+=cart[2] 
# set up initial points on the great circles
    V,XV=[],[]
    if n_planes !=0:
        if n_lines==0:
            V=dir2cart([180.,-45.,1.]) # set the initial direction arbitrarily
        else:
           R=numpy.sqrt(E[0]**2+E[1]**2+E[2]**2) 
           for c in E:
               V.append(c/R) # set initial direction as mean of lines
        U=E[:]   # make a copy of E
        for pole in L:
            XV.append(vclose(pole,V)) # get some points on the great circle
            for c in range(3):
               U[c]=U[c]+XV[-1][c]
# iterate to find best agreement
        angle_tol=1.
        while angle_tol > 0.1:
            angles=[]
            for k in range(n_planes): 
               for c in range(3): U[c]=U[c]-XV[k][c]
               R=numpy.sqrt(U[0]**2+U[1]**2+U[2]**2)
               for c in range(3):V[c]=U[c]/R
               XX=vclose(L[k],V)
               ang=XX[0]*XV[k][0]+XX[1]*XV[k][1]+XX[2]*XV[k][2]
               angles.append(numpy.arccos(ang)*180./numpy.pi)
               for c in range(3):
                   XV[k][c]=XX[c]
                   U[c]=U[c]+XX[c]
               amax =-1
               for ang in angles:
                   if ang > amax:amax=ang
               angle_tol=amax
# calculating overall mean direction and R
        U=E[:]
        for dir in XV:
            for c in range(3):U[c]=U[c]+dir[c]
        R=numpy.sqrt(U[0]**2+U[1]**2+U[2]**2)
        for c in range(3):U[c]=U[c]/R
# get dec and inc of solution points on gt circles
        dirV=cart2dir(U)
# calculate modified Fisher stats fo fit
        n_total=n_lines+n_planes
        NP=n_lines+0.5*n_planes
        if NP<1.1:NP=1.1
        if n_total-R !=0:
            K=(NP-1.)/(n_total-R)
            fac=(20.**(1./(NP-1.))-1.)
            fac=fac*(NP-1.)/K
            a=1.-fac/R
            a95=a
            if abs(a) > 1.0: a95=1.
            if a<0:a95=-a95
            a95=numpy.arccos(a95)*180./numpy.pi
        else: 
            a95=0.
            K='inf'
    else:
        dir=fisher_mean(fdata)
        n_total,R,K,a95=dir["n"],dir["r"],dir["k"],dir["alpha95"]
        dirV[0],dirV[1]=dir["dec"],dir["inc"]
    fpars["tilt_correction"]=tc
    fpars["n_total"]='%i '% (n_total)
    fpars["n_lines"]='%i '% (n_lines)
    fpars["n_planes"]='%i '% (n_planes)
    fpars["R"]='%5.4f '% (R)
    if K!='inf':
        fpars["K"]='%6.0f '% (K)
    else:
        fpars["K"]=K
    fpars["alpha95"]='%7.1f '% (a95)
    fpars["dec"]='%7.1f '% (dirV[0])
    fpars["inc"]='%7.1f '% (dirV[1])
    return fpars

def dobingham(data):
    """
    gets bingham parameters for data
    """
    control,X,bpars=[],[],{}
    N=len(data)
    if N <2:
       return bpars
#
#  get cartesian coordinates
#
    for rec in data:
        X.append(dir2cart([rec[0],rec[1],1.]))
#
#   put in T matrix
#
    T=numpy.array(Tmatrix(X))
    t,V=tauV(T)
    w1,w2,w3=t[2],t[1],t[0]
    k1,k2=binglookup(w1,w2)
    PDir=cart2dir(V[0])
    EDir=cart2dir(V[1])
    ZDir=cart2dir(V[2])
    if PDir[1] < 0: 
        PDir[0]+=180.
        PDir[1]=-PDir[1]
    PDir[0]=PDir[0]%360. 
    bpars["dec"]=PDir[0]
    bpars["inc"]=PDir[1]
    bpars["Edec"]=EDir[0]
    bpars["Einc"]=EDir[1]
    bpars["Zdec"]=ZDir[0]
    bpars["Zinc"]=ZDir[1]
    bpars["n"]=N
#
#  now for Bingham ellipses.
#
    fac1,fac2=-2*N*(k1)*(w3-w1),-2*N*(k2)*(w3-w2)
    sig31,sig32=numpy.sqrt(1./fac1), numpy.sqrt(1./fac2)
    bpars["Zeta"],bpars["Eta"]=2.45*sig31*180./numpy.pi,2.45*sig32*180./numpy.pi
    return  bpars

def doflip(dec,inc):
   """
   flips lower hemisphere data to upper hemisphere
   """
   if inc <0:
       inc=-inc
       dec=(dec+180.)%360.
   return dec,inc

def doincfish(inc):
    """
    gets fisher mean inc from inc only data
    """
    rad,SCOi,SSOi=numpy.pi/180.,0.,0. # some definitions
    abinc=[]
    for i in inc:abinc.append(abs(i))
    MI,std=gausspars(abinc) # get mean inc and standard deviation
    fpars={}
    N=len(inc)  # number of data
    fpars['n']=N
    fpars['ginc']=MI
    if MI<30:
        fpars['inc']=MI
        fpars['k']=0 
        fpars['alpha95']=0 
        fpars['csd']=0 
        fpars['r']=0 
        print 'WARNING: mean inc < 30, returning gaussian mean'
        return fpars
    for i in inc:  # sum over all incs (but take only positive inc)
        coinc=(90.-abs(i))*rad
        SCOi+= numpy.cos(coinc)
        SSOi+= numpy.sin(coinc)
    Oo=(90.0-MI)*rad # first guess at mean
    SCFlag = -1  # sign change flag
    epsilon = float(N)*numpy.cos(Oo) # RHS of zero equations
    epsilon+= (numpy.sin(Oo)**2-numpy.cos(Oo)**2)*SCOi
    epsilon-= 2.*numpy.sin(Oo)*numpy.cos(Oo)*SSOi
    while SCFlag < 0: # loop until cross zero
        if MI > 0 : Oo-=(.01*rad)  # get steeper
        if MI < 0 : Oo+=(.01*rad)  # get shallower
        prev=epsilon
        epsilon = float(N)*numpy.cos(Oo) # RHS of zero equations
        epsilon+= (numpy.sin(Oo)**2.-numpy.cos(Oo)**2.)*SCOi
        epsilon-= 2.*numpy.sin(Oo)*numpy.cos(Oo)*SSOi
        if abs(epsilon) > abs(prev): MI=-1*MI  # reverse direction
        if epsilon*prev < 0: SCFlag = 1 # changed sign
    S,C=0.,0.  # initialize for summation
    for i in inc:
        coinc=(90.-abs(i))*rad
        S+= numpy.sin(Oo-coinc)
        C+= numpy.cos(Oo-coinc)
    k=(N-1.)/(2.*(N-C))
    Imle=90.-(Oo/rad)
    fpars["inc"]=Imle
    fpars["r"],R=2.*C-N,2*C-N
    fpars["k"]=k
    f=fcalc(2,N-1)
    a95= 1. - (0.5)*(S/C)**2 - (f/(2.*C*k))
#    b=20.**(1./(N-1.)) -1.
#    a=1.-b*(N-R)/R
#    a95=numpy.arccos(a)*180./numpy.pi
    csd=81./numpy.sqrt(k)
    fpars["alpha95"]=a95
    fpars["csd"]=csd
    return fpars

def dokent(data,NN):
    """
    gets Kent parameters for data
    """
    X,kpars=[],{}
    N=len(data)
    if N <2:
       return kpars
#
#  get fisher mean and convert to co-inclination (theta)/dec (phi) in radians
#
    fpars=fisher_mean(data)
    pbar=fpars["dec"]*numpy.pi/180.
    tbar=(90.-fpars["inc"])*numpy.pi/180.
#
#   initialize matrices
#
    H=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
    w=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
    b=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
    gam=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
    xg=[]
#
#  set up rotation matrix H
#
    H=[ [numpy.cos(tbar)*numpy.cos(pbar),-numpy.sin(pbar),numpy.sin(tbar)*numpy.cos(pbar)],[numpy.cos(tbar)*numpy.sin(pbar),numpy.cos(pbar),numpy.sin(pbar)*numpy.sin(tbar)],[-numpy.sin(tbar),0.,numpy.cos(tbar)]]
#
#  get cartesian coordinates of data
#
    for rec in data:
        X.append(dir2cart([rec[0],rec[1],1.]))
#
#   put in T matrix
#
    T=Tmatrix(X)
    for i in range(3):
        for j in range(3):
            T[i][j]=T[i][j]/float(N)
#
# compute B=H'TH
#
    for i in range(3):
        for j in range(3):
            for k in range(3):
                w[i][j]+=T[i][k]*H[k][j]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                b[i][j]+=H[k][i]*w[k][j]
#
# choose a rotation w about North pole to diagonalize upper part of B
#
    psi = 0.5*numpy.arctan(2.*b[0][1]/(b[0][0]-b[1][1]))
    w=[[numpy.cos(psi),-numpy.sin(psi),0],[numpy.sin(psi),numpy.cos(psi),0],[0.,0.,1.]]
    for i in range(3):
        for j in range(3):
            gamtmp=0.
            for k in range(3):
                gamtmp+=H[i][k]*w[k][j]      
            gam[i][j]=gamtmp
    for i in range(N):
        xg.append([0.,0.,0.])
        for k in range(3):  
            xgtmp=0.
            for j in range(3):
                xgtmp+=gam[j][k]*X[i][j]
            xg[i][k]=xgtmp
# compute asymptotic ellipse parameters
#
    xmu,sigma1,sigma2=0.,0.,0.
    for  i in range(N):
        xmu+= xg[i][2]
        sigma1=sigma1+xg[i][1]**2
        sigma2=sigma2+xg[i][0]**2
    xmu=xmu/float(N)
    sigma1=sigma1/float(N)
    sigma2=sigma2/float(N)
    g=-2.0*numpy.log(0.05)/(float(NN)*xmu**2)
    if numpy.sqrt(sigma1*g)<1:zeta=numpy.arcsin(numpy.sqrt(sigma1*g))
    if numpy.sqrt(sigma2*g)<1:eta=numpy.arcsin(numpy.sqrt(sigma2*g))
    if numpy.sqrt(sigma1*g)>=1.:zeta=numpy.pi/2.
    if numpy.sqrt(sigma2*g)>=1.:eta=numpy.pi/2.
#
#  convert Kent parameters to directions,angles
#
    kpars["dec"]=fpars["dec"]
    kpars["inc"]=fpars["inc"]
    kpars["n"]=NN
    ZDir=cart2dir([gam[0][1],gam[1][1],gam[2][1]])
    EDir=cart2dir([gam[0][0],gam[1][0],gam[2][0]])
    kpars["Zdec"]=ZDir[0]
    kpars["Zinc"]=ZDir[1]
    kpars["Edec"]=EDir[0]
    kpars["Einc"]=EDir[1]
    if kpars["Zinc"]<0:
        kpars["Zinc"]=-kpars["Zinc"]
        kpars["Zdec"]=(kpars["Zdec"]+180.)%360.
    if kpars["Einc"]<0:
        kpars["Einc"]=-kpars["Einc"]
        kpars["Edec"]=(kpars["Edec"]+180.)%360.
    kpars["Zeta"]=zeta*180./numpy.pi
    kpars["Eta"]=eta*180./numpy.pi
    return kpars

def fshdev(k):
    """
    a call to fshdev(k), where k is kappa, returns a direction from distribution with mean declination of 0, inclination of 90 and kappa of k
    """
    R1=random.random()
    R2=random.random()
    L=numpy.exp(-2*k)
    a=R1*(1-L)+L
    fac=numpy.sqrt((-numpy.log(a))/(2*k))
    inc=90.-2*numpy.arcsin(fac)*180./numpy.pi
    dec=2*numpy.pi*R2*180./numpy.pi
    return dec,inc

def fishrot(kappa,N,D,I):
    """
    Description: generates set of Fisher distributed data from specified distribution 
	Input: kappa (fisher distribution concentration parameter), number of desired subsamples, Dec and Inc
	Output: list with N pairs of Dec, Inc.
    """
    out_d=[]
    out=[]

    for k in range(N): 
        dec,inc= fshdev(kappa)  # send kappa to fshdev
        drot,irot=dodirot(dec,inc,D,I)
        out_d=[drot,irot]
        out.append(out_d)
    return out
	
def dimap(D,I):
    """
    FUNCTION TO MAP DECLINATION, INCLINATIONS INTO EQUAL AREA PROJECTION, X,Y

    Usage:     dimap(D, I)
    Argin:     Declination (float) and Inclination (float)

    """
### DEFINE FUNCTION VARIABLES
    XY=[0.,0.]                                     # initialize equal area projection x,y

### GET CARTESIAN COMPONENTS OF INPUT DIRECTION
    X=dir2cart([D,I,1.])

### CHECK IF Z = 1 AND ABORT
    if X[2] ==1.0: return XY                       # return [0,0]

### TAKE THE ABSOLUTE VALUE OF Z
    if X[2]<0:X[2]=-X[2]                           # this only works on lower hemisphere projections

### CALCULATE THE X,Y COORDINATES FOR THE EQUAL AREA PROJECTION
    R=numpy.sqrt( 1.-X[2])/(numpy.sqrt(X[0]**2+X[1]**2)) # from Collinson 1983
    XY[1],XY[0]=X[0]*R,X[1]*R

### RETURN XY[X,Y]
    return XY

def dimap_V(D,I):
    """
    FUNCTION TO MAP DECLINATION, INCLINATIONS INTO EQUAL AREA PROJECTION, X,Y

    Usage:     dimap_V(D, I)
        D and I are both numpy arrays

    """
### GET CARTESIAN COMPONENTS OF INPUT DIRECTION
    DI=numpy.array([D,I]).transpose() # 
    X=dir2cart(DI).transpose()
### CALCULATE THE X,Y COORDINATES FOR THE EQUAL AREA PROJECTION
    R=numpy.sqrt( 1.-abs(X[2]))/(numpy.sqrt(X[0]**2+X[1]**2)) # from Collinson 1983
    XY=numpy.array([X[1]*R,X[0]*R]).transpose()

### RETURN XY[X,Y]
    return XY
	
def gaussdev(mean,sigma):
    """
    returns a number randomly drawn from a gaussian distribution with the given mean, sigma
    """
    return random.normal(mean,sigma) # return gaussian deviate

def get_unf(N):
    """
    Called with get_unf(N).
 subroutine to retrieve N uniformly distributed directions
 using the way described in Fisher et al. (1987).
    """
#
# get uniform directions  [dec,inc]
    z=random.uniform(-1.,1.,size=N)
    t=random.uniform(0.,360.,size=N) # decs
    i=numpy.arcsin(z)*180./numpy.pi # incs
    return numpy.array([t,i]).transpose()

#def get_unf(N): #Jeff's way
    """
     subroutine to retrieve N uniformly distributed directions
    """

def cross(v, w):
    """
     cross product of two vectors
    """
    x = v[1]*w[2] - v[2]*w[1]
    y = v[2]*w[0] - v[0]*w[2]
    z = v[0]*w[1] - v[1]*w[0]
    return [x, y, z]

def apseudo(Ss,ipar,sigma):
    """
     draw a bootstrap sample of Ss
    """
#
    Is=random.randint(0,len(Ss)-1,size=len(Ss)) # draw N random integers
    Ss=numpy.array(Ss)
    if ipar==0:
        BSs=Ss[Is]
    else: # need to recreate measurement - then do the parametric stuffr
        A,B=design(6) # get the design matrix for 6 measurements
        K,BSs=[],[]
        for k in range(len(Ss)):
            K.append(numpy.dot(A,Ss[k]))
        Pars=numpy.random.normal(K,sigma)
        for k in range(len(Ss)):
            BSs.append(numpy.dot(B,Pars[k]))
    return numpy.array(BSs)

def sbootpars(Taus,Vs):
    """
     get bootstrap parameters for s data
    """
#
    Tau1s,Tau2s,Tau3s=[],[],[]
    V1s,V2s,V3s=[],[],[]
    nb=len(Taus)
    bpars={}
    for k in range(nb):
        Tau1s.append(Taus[k][0])
        Tau2s.append(Taus[k][1])
        Tau3s.append(Taus[k][2])
        V1s.append(Vs[k][0])
        V2s.append(Vs[k][1])
        V3s.append(Vs[k][2])
    x,sig=gausspars(Tau1s) 
    bpars["t1_sigma"]=sig
    x,sig=gausspars(Tau2s) 
    bpars["t2_sigma"]=sig
    x,sig=gausspars(Tau3s) 
    bpars["t3_sigma"]=sig
    kpars=dokent(V1s,len(V1s))
    bpars["v1_dec"]=kpars["dec"]
    bpars["v1_inc"]=kpars["inc"]
    bpars["v1_zeta"]=kpars["Zeta"]*numpy.sqrt(nb)
    bpars["v1_eta"]=kpars["Eta"]*numpy.sqrt(nb)
    bpars["v1_zeta_dec"]=kpars["Zdec"]
    bpars["v1_zeta_inc"]=kpars["Zinc"]
    bpars["v1_eta_dec"]=kpars["Edec"]
    bpars["v1_eta_inc"]=kpars["Einc"]
    kpars=dokent(V2s,len(V2s))
    bpars["v2_dec"]=kpars["dec"]
    bpars["v2_inc"]=kpars["inc"]
    bpars["v2_zeta"]=kpars["Zeta"]*numpy.sqrt(nb)
    bpars["v2_eta"]=kpars["Eta"]*numpy.sqrt(nb)
    bpars["v2_zeta_dec"]=kpars["Zdec"]
    bpars["v2_zeta_inc"]=kpars["Zinc"]
    bpars["v2_eta_dec"]=kpars["Edec"]
    bpars["v2_eta_inc"]=kpars["Einc"]
    kpars=dokent(V3s,len(V3s))
    bpars["v3_dec"]=kpars["dec"]
    bpars["v3_inc"]=kpars["inc"]
    bpars["v3_zeta"]=kpars["Zeta"]*numpy.sqrt(nb)
    bpars["v3_eta"]=kpars["Eta"]*numpy.sqrt(nb)
    bpars["v3_zeta_dec"]=kpars["Zdec"]
    bpars["v3_zeta_inc"]=kpars["Zinc"]
    bpars["v3_eta_dec"]=kpars["Edec"]
    bpars["v3_eta_inc"]=kpars["Einc"]
    return bpars

def s_boot(Ss,ipar,nb):
    """
     returns bootstrap parameters for S data
    """
    npts=len(Ss)
# get average s for whole dataset
    nf,Sigma,avs=sbar(Ss)
    Tmean,Vmean=doseigs(avs) # get eigenvectors of mean tensor
#
# now do bootstrap to collect Vs and taus of bootstrap means
#
    Taus,Vs=[],[]  # number of bootstraps, list of bootstrap taus and eigenvectors
#

    for k in range(nb): # repeat nb times
#        if k%50==0:print k,' out of ',nb
        BSs=apseudo(Ss,ipar,Sigma) # get a pseudosample - if ipar=1, do a parametric bootstrap
        nf,sigma,avbs=sbar(BSs) # get bootstrap mean s
        tau,Vdirs=doseigs(avbs) # get bootstrap eigenparameters
        Taus.append(tau)
        Vs.append(Vdirs)
    return Tmean,Vmean,Taus,Vs

def pseudo(DIs):
    """
     draw a bootstrap sample of Directions
    """
#
    Inds=numpy.random.randint(len(DIs),size=len(DIs))
    D=numpy.array(DIs)
    return D[Inds]

def di_boot(DIs):
    """
     returns bootstrap parameters for Directional data
    """
# get average DI for whole dataset
    fpars=fisher_mean(DIs)
#
# now do bootstrap to collect BDIs  bootstrap means
#
    nb,BDIs=5000,[]  # number of bootstraps, list of bootstrap directions
#
    
    for k in range(nb): # repeat nb times
#        if k%50==0:print k,' out of ',nb
        pDIs= pseudo(DIs) # get a pseudosample 
        bfpars=fisher_mean(pDIs) # get bootstrap mean bootstrap sample
        BDIs.append([bfpars['dec'],bfpars['inc']])
    return BDIs

def pseudosample(x):
    """
     draw a bootstrap sample of x
    """
#
    BXs=[]
    for k in range(len(x)):
        ind=random.randint(0,len(x)-1)
        BXs.append(x[ind])
    return BXs 
