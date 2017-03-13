print '\nThis program uses the PmagPy software draw utilities \n(Tauxe et al. 2016, G3, DOI: 10.1002/2016GC006307)'

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv


rad=np.pi/180
deg=180/np.pi

def eqareaXY(dato):
    d=dato[0]
    i=dato[1]
    x=(2.*np.sin(rad*(90.-np.absolute(i))/2.)*np.sin(rad*d))/np.sqrt(2)
    y=(2.*np.sin(rad*(90.-np.absolute(i))/2.)*np.cos(rad*d))/np.sqrt(2)
    XY=[x,y]
    
    return x,y

def plot_net(fignum):#From Tauxe et al(2016)
    """
    Draws circle and tick marks for equal area projection.
    """

# make the perimeter
    plt.figure(num=fignum)
    plt.clf()
    plt.axis("off")
    Dcirc=np.arange(0,361.)
    Icirc=np.zeros(361,'f')
    Xcirc,Ycirc=[],[]
    for k in range(361):
        XY=dimap(Dcirc[k],Icirc[k])
        Xcirc.append(XY[0])
        Ycirc.append(XY[1])
    plt.plot(Xcirc,Ycirc,'k')

# put on the tick marks
    Xsym,Ysym=[],[]
    for I in range(10,100,10):
        XY=dimap(0.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    plt.plot(Xsym,Ysym,'k+')
    Xsym,Ysym=[],[]
    for I in range(10,90,10):
        XY=dimap(90.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    plt.plot(Xsym,Ysym,'k+')
    Xsym,Ysym=[],[]
    for I in range(10,90,10):
        XY=dimap(180.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    plt.plot(Xsym,Ysym,'k+')
    Xsym,Ysym=[],[]
    for I in range(10,90,10):
        XY=dimap(270.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    plt.plot(Xsym,Ysym,'k+')
    for D in range(0,360,10):
        Xtick,Ytick=[],[]
        for I in range(4):
            XY=dimap(D,I)
            Xtick.append(XY[0])
            Ytick.append(XY[1])
        plt.plot(Xtick,Ytick,'k')
    plt.axis("equal")
    plt.axis((-1.05,1.05,-1.05,1.05))
    
def dimap(D,I):#From Tauxe et al(2016)
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
    R=np.sqrt( 1.-X[2])/(np.sqrt(X[0]**2+X[1]**2)) # from Collinson 1983
    XY[1],XY[0]=X[0]*R,X[1]*R

### RETURN XY[X,Y]
    return XY

def circ(dec,dip,alpha):#From Tauxe et al(2016)
    """
    function to calculate points on an circle about dec,dip with angle alpha
    """
    rad=np.pi/180.
    D_out,I_out=[],[]
    dec,dip,alpha=dec*rad ,dip*rad,alpha*rad
    dec1=dec+np.pi/2.
    isign=1
    if dip!=0: isign=(abs(dip)/dip)
    dip1=(dip-isign*(np.pi/2.))
    t=[[0,0,0],[0,0,0],[0,0,0]]
    v=[0,0,0]
    t[0][2]=np.cos(dec)*np.cos(dip)
    t[1][2]=np.sin(dec)*np.cos(dip)
    t[2][2]=np.sin(dip)
    t[0][1]=np.cos(dec)*np.cos(dip1)
    t[1][1]=np.sin(dec)*np.cos(dip1)
    t[2][1]=np.sin(dip1)
    t[0][0]=np.cos(dec1)
    t[1][0]=np.sin(dec1)
    t[2][0]=0
    for i in range(101):
        psi=float(i)*np.pi/50.
        v[0]=np.sin(alpha)*np.cos(psi)
        v[1]=np.sin(alpha)*np.sin(psi)
        v[2]=np.sqrt(abs(1.-v[0]**2 - v[1]**2))
        elli=[0,0,0]
        for j in range(3):
            for k in range(3):
                elli[j]=elli[j] + t[j][k]*v[k]
        Dir=cart2dir(elli)
        D_out.append(Dir[0])
        I_out.append(Dir[1])
    return D_out,I_out
    
def dir2cart(d):
   # converts list or array of vector directions, in degrees, to array of cartesian coordinates, in x,y,z
    ints=np.ones(len(d)).transpose() # get an array of ones to plug into dec,inc pairs
    d=np.array(d)
    rad=np.pi/180.
    if len(d.shape)>1: # array of vectors
        decs,incs=d[:,0]*rad,d[:,1]*rad
        if d.shape[1]==3: ints=d[:,2] # take the given lengths
    else: # single vector
        decs,incs=np.array(d[0])*rad,np.array(d[1])*rad
        if len(d)==3:
            ints=np.array(d[2])
        else:
            ints=np.array([1.])
    cart= np.array([ints*np.cos(decs)*np.cos(incs),ints*np.sin(decs)*np.cos(incs),ints*np.sin(incs)]).transpose()
    return cart
 
def cart2dir(cart):
    """
    converts a direction to cartesian coordinates.  takes an array of [x,y,z])
    """
    cart=np.array(cart)
    rad=np.pi/180. # constant to convert degrees to radians
    if len(cart.shape)>1:
        Xs,Ys,Zs=cart[:,0],cart[:,1],cart[:,2]
    else: #single vector
        Xs,Ys,Zs=cart[0],cart[1],cart[2]
    Rs=np.sqrt(Xs**2+Ys**2+Zs**2) # calculate resultant vector length
    Decs=(np.arctan2(Ys,Xs)/rad)%360. # calculate declination taking care of correct quadrants (arctan2) and making modulo 360.
    try:
        Incs=np.arcsin(Zs/Rs)/rad # calculate inclination (converting to degrees) #
    except:
        print 'trouble in cart2dir' # most likely division by zero somewhere
        return np.zeros(3)

    return np.array([Decs,Incs,Rs]).transpose() # return the directions list
 
def saveInputFile(files):
    '''
    Get information coming from a *txt file (files) in a list of list without its header.
    Input: string with the name of the file
           file: (Site, Dec, Inc, alfa95, dip direction, strike, dip)
    Output: (i) data: list of list with Site, Dec, Inc, alfa95, strike, dip (without header)
            (ii) geo: list of list with values of magnetization GEO; Dec, Inc, alfa95 (without header)
            (iii) bed: list of list with bedding data; Strike (RHR) and Dip
    
    '''
    
    #beig 'files' a csv comma separated with 5 columns and header (Site, Dec, Inc, alfa95, bed strike)   #Return a list with the same but without header, and with the numbers as float
    #Saving the input file/data 
    #output have the same colums, but whithou header.
    reader=csv.reader(open(files, 'rU'), delimiter=' ')
    dat=list(reader)
    file.close()
    data=dat[1:]    #removing the header
    data_float=data[:]

    #Converting the input data to float
    for i, lista in enumerate(data):
        for j, val in enumerate(lista):
            if j>0:   #condicional es para dejar el nombre del site tranquilo
                data_float[i][j]=float(val)
                
    geo=[]
    tilt=[]
    site=[]
    sc=[]
    bfd=[]

    for dato in data_float:
        site_site=[dato[0]]
        site.append(site_site)
        dif=(dato[10]-dato[1])%360
        if dif>90 and dif<270: 
            st=(dato[10]+180)%360
        else: st=dato[10]
        sc_site=[st,0.,dato[11]]
        sc.append(sc_site)
        geo_site=[dato[1],dato[2], dato[3]]
        geo.append(geo_site)
        tilt_site=[dato[6],dato[7], dato[3]]
        tilt.append(tilt_site)
        bfd_site=[dato[8],dato[9], dato[3]]
        bfd.append(bfd_site)
    
    return site,sc,geo,tilt,bfd
 
def smallcirc(a,zorder=1):#From Tauxe et al(2016)
    Ds,Is=circ(a[0],a[1],a[2])
    Xcirc,Ycirc=[],[]
    for k in range(len(Ds)):
        XY=dimap(Ds[k],Is[k])
        Xcirc.append(XY[0])
        Ycirc.append(XY[1])
    plt.plot(Xcirc,Ycirc,'k',zorder=zorder)

def plot_di_mean(dec,inc,a95,color='k',marker='o',markersize=20,label='',legend='no',zorder=2):#From Tauxe et al(2016)
    """
    Plot a mean direction (declination, inclination) with alpha_95 ellipse on
    an equal area plot.

    Before this function is called, a plot needs to be initialized with code
    that looks something like:
    >fignum = 1
    >plt.figure(num=fignum,figsize=(10,10),dpi=160)
    >ipmag.plot_net(fignum)

    Required Arguments
    -----------
    dec : declination of mean being plotted
    inc : inclination of mean being plotted
    a95 : a95 confidence ellipse of mean being plotted

    Optional Keywords
    -----------
    color : the default color is black. Other colors can be chosen (e.g. 'r').
    marker : the default is a circle. Other symbols can be chosen (e.g. 's').
    markersize : the default is 20. Other sizes can be chosen.
    label : the default is no label. Labels can be assigned.
    legend : the default is no legend ('no'). Putting 'yes' will plot a legend.
    """
    DI_dimap=dimap(dec,inc)
    if inc < 0:
        plt.scatter(DI_dimap[0],DI_dimap[1],
        edgecolors=color ,facecolors='white',
        marker=marker,s=markersize,label=label)
    if inc >= 0:
        plt.scatter(DI_dimap[0],DI_dimap[1],
        edgecolors=color,facecolors=color,
        marker=marker,s=markersize,label=label)
    Xcirc,Ycirc=[],[]
    Da95,Ia95=circ(dec,inc,a95)
    if legend=='yes':
        plt.legend(loc=2)
    for k in range(len(Da95)):
        XY=dimap(Da95[k],Ia95[k])
        Xcirc.append(XY[0])
        Ycirc.append(XY[1])
    plt.plot(Xcirc,Ycirc,c=color,zorder=zorder)
    plt.tight_layout()

print '\nFile input is the main file output from the pySCIs_calc.py program'
files = raw_input("File name: ")
try:
    file = open(files)
except:
    print "Try again (do not forget the extension of the file)"
    exit() 
 
while True:
    try:
        DecInp = float(raw_input("Declination of the expected direction: "))
        if DecInp <0 or DecInp>360:
            print 'Declination should be in range 0-360'
            sys.exit(1)
        break
    except ValueError:
        print "This is not a valid declination"

while True:
    try:
        IncInp = float(raw_input("Inclination of the expected direction: "))
        #if DecInp <-90 or DecInp>90:
        #    print 'Declination should be in range +-90'
        #    sys.exit(1)
        break
    except ValueError:
        print "This is not a valid inclination"
exp_dir=[DecInp, IncInp] 
 
site,sc,geo,tilt,bfd=saveInputFile(files)   #Saving in diferent list the data
name_svg=files[:-9]+'_labels.svg'
name_png=files[:-9]+'_labels.png'
n=len(site)
 
plt.figure(num=1,figsize=(8,8))
plot_net(fignum=1)

plt.text(0.76, 0.97, 'n='+str(n), fontsize = 20)
plt.text(0.8, 0.9, 'BBC', fontsize = 15)
plt.text(0.8, 0.8, 'ATBC', fontsize = 15)
plt.text(0.8, 0.7, 'BFD', fontsize = 15)
plt.scatter(0.75, 0.91, color='r',marker='s',s=40)
plt.scatter(0.75, 0.81, color='g',marker='^',s=60)
plt.scatter(0.75, 0.71, color='b',marker='o',s=40)
plt.text(-0.87, 0.92, 'Reference', fontsize = 15)
plt.scatter(-0.95, 0.95, color='m',marker='*',s=1000,)

plot_di_mean(exp_dir[0],exp_dir[1],0,color='m',marker='*',markersize=1000,label='Reference',legend='no',zorder=4)

for dato in sc:
    smallcirc(dato,zorder=1)

for dato in geo:
    plot_di_mean(dato[0],dato[1],dato[2],color='r',marker='s',markersize=40,label='Geo',legend='no',zorder=3)

for dato in tilt:
    plot_di_mean(dato[0],dato[1],dato[2],color='g',marker='^',markersize=50,label='Tilt',legend='no',zorder=3)

for dato in bfd:
    plot_di_mean(dato[0],dato[1],dato[2],color='b',marker='o',markersize=40,label='BFD',legend='no',zorder=3)
    
for dato_t, dato_s in zip (tilt,site):
    x,y=eqareaXY(dato_t)
    plt.text(x+0.03, y-0.03, dato_s[0], fontsize = 10, zorder=4)

plt.savefig(name_svg,transparent='True')
plt.savefig(name_png,transparent='True',dpi=300)
print '\nThe file was saved as', name_svg, 'and as', name_png
plt.show()
quit()
