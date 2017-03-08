# -*- coding: utf-8 -*-

# In[1]:

import csv
import numpy as np
import sys
rad=np.pi/180
deg=180/np.pi
from time import time


# In[2]:

def saveInputFile(files):
    '''
    Save a *csv file (files) in a list of list without its header.
    Input: string with the name of the file
           file: (Site, Dec, Inc, alfa95, dip direction, strike, dip)
    Output: (i) data: list of list with Site, Dec, Inc, alfa95, strike, dip (without header)
            (ii) geo: list of list with values of magnetization GEO; Dec, Inc, alfa95 (without header)
            (iii) bed: list of list with bedding data; Strike (RHR) and Dip
    
    '''
    
    #beig 'files' a csv comma separated with 5 columns and header (Site, Dec, Inc, alfa95, bed strike)   #Return a list with the same but without header, and with the numbers as float
    #Saving the input file/data 
    #output have the same colums, but whithou header.
    #file = open(files)
    reader=csv.reader(open(files, 'rU'), delimiter=' ')
    dat=list(reader)
    file.close()
    data=dat[1:]    #guardamos los datos sin el encabezado
    data_float=data[:]

    #Converting the input data to float
    #Es posible que esta conversion de problemas. En ese caso, hay un modulo que se llama decimal que se puede utilizar, o quiza decirle que lo redondee. Aunque creo que lo hacia bien y no ponia los decimales que le daba la gana, sino los que traian los datos.
    for i, lista in enumerate(data):
        for j, val in enumerate(lista):
            if j>0:   #condicional es para dejar el nombre del site tranquilo
                data_float[i][j]=float(val)
    
    geo=[]
    bed=[]
    site=[]
    #print data_float

    for dato in data_float:
        site_site=[dato[0]]
        site.append(site_site)
        geo_site=[dato[1],dato[2], dato[3]]
        geo.append(geo_site)
        bed_site=[((dato[4]-90)%360),dato[5]]
        bed.append(bed_site)
    
    return data_float, site, geo, bed


# In[3]:

def dir2cart(dir): #being 'dir' a list of list, with pairs of Dec,Inc in degrees
    #changes from geographic to cartesians coord
    cart=[]   #output is a list of list, with x,y,z
    rad=np.pi/180

    for site in dir:
        x=np.cos(rad*site[0])*np.cos(rad*site[1]) 
        y=np.sin(rad*site[0])*np.cos(rad*site[1]) 
        z=np.sin(rad*site[1])                 
        cartSite=[x,y,z]
        cart.append(cartSite)
    return cart


# In[4]:

def cart2dir(cart): #being 'cart' a list of list, with x, y, z
    ##changes from geographic to cartesians coord
    dir=[]    #output is a list of list, with Dec and Inc (in degree)
    deg=180/np.pi
    
    for site in cart:
        Dec=round(deg*np.arctan2(site[1],site[0]),1)
        Inc=round(deg*np.arcsin(site[2]), 1)
        dirSite=[Dec,Inc]
        dir.append(dirSite)
    return dir


# In[5]:

def fisher_mean(dir):
    '''
    It does fisher mean
        Input: list of list, with pairs of Dec, Inc (degrees)
        Output: a list with [Dec_mean, Inc_mean, kappa, alfa95]
    '''
    
    fisher_mean=[]
    
    #change to cartesian coordinats using dir2cart defined funcion
    car=dir2cart(dir)
    N=len(dir)
    X,Y,Z=0,0,0
    
    for site in car:
        X+=site[0]
        Y+=site[1]
        Z+=site[2]

    R=np.sqrt(X**2+Y**2+Z**2)
    k=(N-1.)/(N-R)
    
    cosAlfa=1.-((N-R)/R)*(20.**(1./(N-1.))-1.)
    if cosAlfa<-1:cosAlfa=-1
    alpha95=deg*np.arccos(cosAlfa)

    #calculating mean direction
    DecMean=deg*np.arctan(Y/X)
    if X<0.:
        DecMean=DecMean+180.
    elif Y<0.:
        DecMean=DecMean+360.
    IncMean=deg*np.arcsin(Z/R)
    fisher_mean=[DecMean,IncMean,k,alpha95]
    return fisher_mean


# In[6]:

def tilt_rot(m,b):
    """
    This function allow to apply the bedding correction to a vector
    Input data: list of list with [m] the magnetic vector ( Dec and Inc in degrees) 
                                  [b] the corresponding bed strike (RHR) and the dip of the bed
                                  len(m)=len(b)
    Ouput data: list of list with the rotated vectors
    """
        
    rotated=[]
    
    
    for datom, datob in zip(m,b):
        D=rad*datom[0]   #declination
        I=rad*datom[1]    #inclination
        s=rad*datob[0]    #strike
        d=rad*datob[1]    #dip
        dif=D-s
                  
        x=np.sin(s)*np.cos(dif)*np.cos(I)+np.cos(s)*np.cos(d)*np.sin(dif)*np.cos(I)+np.cos(s)*np.sin(d)*np.sin(I)
        y=np.cos(s)*np.cos(dif)*np.cos(I)-np.sin(s)*np.cos(d)*np.sin(dif)*np.cos(I)-np.sin(s)*np.sin(d)*np.sin(I)
        z=-np.sin(d)*np.sin(dif)*np.cos(I)+np.cos(d)*np.sin(I)
        
        Dec=deg*np.arctan(x/y)
        Inc=deg*np.arcsin(z)
        
        if y<0:
            Dec=Dec+180.
        elif x<0:
            Dec=Dec+360.
        
        outSite=[Dec,Inc]
        rotated.append(outSite)

    return rotated


# In[7]:

def paleo_dip(tilt,bed,q):
    """
    This function allow to calculate the paledip
    For this, it rotates 90º the TITL magnetization and the BFD with a axis perpendicular to the strike.
    In this way, the strike will be upright, the unfolding angle is the difference between TILT and BFD and the sense of rotation is clear
    
    Input data: list of list with [tilt] Dec and Inc (in degrees) of the tilt magnetic vector
                                  [bed] strike (RHR) and dip of bedding
                                  [q] Dec and Inc (in degrees) of the BFD for each site
    Ouput data: list of list with the paleobuz (Strike and dip)
    """
        
    paleobed=[]    
    rot=[]
    
    #calculating the rotation axis that allow rotate the strike to the upright position (looking down)
    for dato in bed:
        s=(dato[0]+90.)%360.   #declination
        sit=[s,90.]
        rot.append(sit)
    
    #rotating both tilt and q data.
    tilt_r=tilt_rot(tilt,rot)
    q_r=tilt_rot(q,rot)
    
    #calculation BFD for each site
    for dato_tilt, dato_q, dato_bed in zip(tilt_r,q_r,bed):
        paleo_dip=(dato_q[0]-dato_tilt[0])%360.
        if paleo_dip>=180:
            paleo_st=(dato_bed[0]-180.)%360.
            paleo_dip=360.-paleo_dip
        else:
            paleo_st=dato_bed[0]
        
        sit=[paleo_st,paleo_dip]
        paleobed.append(sit)  
        
    return paleobed


# In[8]:

def calAQQ2(geo,bed,point):
    """
    This function calculates (i)the directions over each small circle closest to a one point
                             (ii)the sum of all minimum angles between one point (P) and all SCs
    Data: is a list of list (without headets) with Site, Dec, Inc, alfa and bed strike (in degrees)
    Output: (i) Ai, the sum of the angles.
            (ii)out, a list of list with the name of the site and all closest directions.
    """
        
    q_rot=[]
    rot=[]
    unrot=[]
    point_l=[]
    A=0
    
    for dato in bed:
        s=(rad*dato[0]-rad*90.)%(2.*np.pi)   #rotation parameters
        sit=[deg*s,90.]
        unsit=[deg*s,-90.]
        rot.append(sit)
        unrot.append(unsit)
        point_l.append(point)
    
    geo_rot=tilt_rot(geo,rot)
    p_rot=tilt_rot(point_l,rot)
    
    for dato_geo, dato_p in zip(geo_rot,p_rot):       
        Dqi=dato_p[0]
        Iqi=dato_geo[1]
        
        alfai=abs(dato_geo[1]-dato_p[1])      
        A += alfai
        
        site=[Dqi,Iqi]
        q_rot.append(site)
        
    q=tilt_rot(q_rot,unrot)
        
    return q,A


# In[9]:

def calA(geo,bed,point):
    """
    This function calculates (i)the sum of all minimum angles between one point (P) and all SCs
    Input: geo: list of list with Dec and Inc (in degrees)

           bed: list of list with Strike and dip
           point: list with the Dec and Inc of one direction
    Output: (i) A, the sum of the angles.

    """
        
    q_rot=[]
    rot=[]
    point_l=[]
    A=0
    
    for dato in bed:
        s=(rad*dato[0]-rad*90.)%(2.*np.pi)   #rotation parameters
        sit=[deg*s,90.]
        rot.append(sit)
        point_l.append(point)
    
    geo_rot=tilt_rot(geo,rot)
    p_rot=tilt_rot(point_l,rot)
    
    for dato_geo, dato_p in zip(geo_rot,p_rot):       
        alfai=abs(dato_geo[1]-dato_p[1])
        A += alfai
        
    return A


# In[10]:

def cal_matriz(geo,bed):
	
	
    matriz_pos=[]
    d_mat=0.
    n=len(geo)
	


    while d_mat<360:
        i_mat=0.000001
        while i_mat<90:
            point=[d_mat,i_mat]
            x=2.*np.sin(rad*(90.-np.absolute(i_mat))/2.)*np.sin(rad*d_mat)
            xx=x/np.sqrt(2)
            y=2.*np.sin(rad*(90.-np.absolute(i_mat))/2.)*np.cos(rad*d_mat)
            yy=y/np.sqrt(2)
            A=round(calA(geo,bed,point),3)
            An=round(A/n,3)
            out=[d_mat,round(i_mat,0),x,y,xx,yy,A,An]
            matriz_pos.append(out)
            i_mat+=1
        d_mat+=1
    d_mat=0
    
    return matriz_pos


# In[11]:

def save_out_file(header,data,name):
    name_csv=name+'.dat'
    data_out=[header]
    data_out.extend(data)
    file_out=open(name_csv, 'w')
    writer=csv.writer(file_out, delimiter=' ', lineterminator='\r')
    for row in data_out:
        writer.writerow(row)
    del writer
    file_out.close()


# In[12]:

def ang2point(a,b):
    """
    It calculates the angle between two points along a great circle.
    The coordinates of the points are introduced as two list; Dec and Inc (degrees).
    Te output is a two decimal float with the angle in degrees
    """

    xa=np.cos(rad*a[0])*np.cos(rad*a[1]) 
    ya=np.sin(rad*a[0])*np.cos(rad*a[1]) 
    za=np.sin(rad*a[1])                 
    a=[xa,ya,za]
    
    xb=np.cos(rad*b[0])*np.cos(rad*b[1]) 
    yb=np.sin(rad*b[0])*np.cos(rad*b[1]) 
    zb=np.sin(rad*b[1])
    b=[xb,yb,zb]
    
    c=deg*np.arccos(np.dot(a,b))
    c=round(c,2)
    return c


# In[49]:

def cal_api(geo,bed):
    """
    The apical angle of the small circle is calculated for each site.
    input: 
        geo: list of list with Dec and Inc (in degrees)
        strike: list of list with strike an dip (only strike is used)
    """
    
    api=[]
    
    for dato_geo, dato_bed in zip(geo,bed):
        api_site=deg*np.arccos(np.cos(rad*(dato_geo[0]-dato_bed[0]))*np.cos(rad*dato_geo[1]))
        api.append(api_site)
        
    return api

def cal_api2(geo,bed):
	"""
	The apical angle of the small circle is calculated for each site.
	input: 
		geo: list of list with Dec and Inc (in degrees)
		strike: list of list with strike an dip (only strike is used)
	output:
		At: a list with the apical angle and the strike (in degrees)
	"""
	
	At=[]
	for dato_geo, dato_bed in zip(geo,bed):
		api_site=deg*np.arccos(np.cos(rad*(dato_geo[0]-dato_bed[0]))*np.cos(rad*dato_geo[1]))
		trend=dato_bed[0]
		at_site=[api_site,trend]
		At.append(at_site)
	return At
	
def intersec(sc1,sc2):
    """
    This function calculates the posible intersection between two give SCs
    input:
        SC1, SC2, a couple of SC with [apical angle, trend]
    output:
        [Dec, Inc]
    """
    d1=np.cos(rad*sc1[0])
    d2=np.cos(rad*sc2[0])
    rt1=rad*sc1[1]
    rt2=rad*sc2[1]
    a=d1*np.cos(rt2)-d2*np.cos(rt1)
    b=d1*np.sin(rt2)-d2*np.sin(rt1)
    D=np.arctan(-a/b)
    I=deg*np.arccos(d1/np.cos(D-rt1))
    D=(deg*D+360.)%360.
    out=[round(D,2),round(I,2)]
    
    return out

def inter(api2,site):
    """
    This function calculates the posible intersection between the SCs
    First, it test if two SCs intersect, and if true, intersec(SC1,SC2) calculates de intersection
    input:
        api2: a list of list with the apical angle and the trend (in degrees) of each site
        
    """
    i=0
    j=1
    l=len(api2)
    
    in_out=[]

    while i<l:
        while j<l:
            difSt=deg*np.arccos(np.cos(rad*(api2[i][1]-api2[j][1]))) #calculando la diferencia de strike
            sumAp=api2[i][0]%90+api2[j][0]%90
            resAP=np.abs(api2[i][0]-api2[j][0])
            if sumAp>=difSt:
                if resAP<=difSt:
                    inter_ij=intersec(api2[i],api2[j])
                    inter_ij.extend(site[i])
                    inter_ij.extend(site[j])
                    in_out.append(inter_ij)
            j+=1
        i+=1
        j=i+1
    return in_out
	
	
def minA(geo,bed,Q,A):
    '''
    This function found the direction which minimize A value
    '''
    angle=2. #inicialiting angle
    Qmean=fisher_mean(Q)
    while angle>0.05:
        Qnew,A=calAQQ2(geo,bed,Qmean)
        Qnew_mean=fisher_mean(Qnew)
        angle=ang2point(Qnew_mean,Qmean)
        Qmean=Qnew_mean[:]
        
    return Qnew_mean,Qnew,A

#The input file
files = raw_input("\nFile name: ")
try:
    file = open(files)
except:
    print "Try again (do not forget the extension of the file)"
    exit()


data,site,geo,bed=saveInputFile(files)   #Saving in diferent list the data
tilt=tilt_rot(geo,bed)                   #Calculating TILT directions

#Interacting... The workflow
pregunta_remag = raw_input("Do you want to calculate the remagnetization direction? (y/n) ")
if pregunta_remag=='y':
    pregunta_matriz= raw_input("Do you want to calculate the A/N matrix? (y/n) ")
else: pregunta_matriz='n'

#Control point of the raw_input
if pregunta_remag<>'y' and pregunta_remag<>'n':
	print "\nTake care, you have to answer 'y' or 'n'."
	sys.exit(1)
if pregunta_matriz<>'y' and pregunta_matriz<>'n':
	print "\nTake care, you have to answer 'y' or 'n'."
	sys.exit(1)

if pregunta_remag=='y':
    pregunta_inc= raw_input("\nAre you looking for a (p)ositive or a (n)egative inclination \nof the remagnetization direction? (p/n) ")
    if pregunta_inc<>'p' and pregunta_inc<>'n':
        print "\nTake care, you have to answer 'p' or 'n'."  
        sys.exit(1)

#The user have to input the remagnetization direction
while True:
    if pregunta_remag=='n':
        DecInp = float(raw_input("Declination of the remagnetization direction: "))
        try:
            if DecInp <0 or DecInp>360:
                print 'Declination should be in range 0-360'
                sys.exit(1)
            break
        except ValueError:
            print "This is not a valid declination"
    break

while True:
    if pregunta_remag=='n':
        IncInp = float(raw_input("Inclination of the remagnetization direction: "))
        try:
            if IncInp <-90. or IncInp>90.:
                print 'Inclination should be in range +-90'
                sys.exit(1)
            break
        except ValueError:
            print "This is not a valid inclination"
    break


#If the remagnetization direction is calculating, four starting points are inicialited	
if pregunta_remag=='n': 
    point=[DecInp, IncInp]
else: 
    if pregunta_inc=='p':
        point1=[60.,20.]
        point2=[180.,20.]
        point3=[300.,20.]
        point4=[60.,90.]
    else:
        point1=[60.,-20.]
        point2=[180.,-20.]
        point3=[300.,-20]
        point4=[60.,-90.]


time_ini=time()


#Calculating the mean of the calculated Q directions.
if pregunta_remag=='n':
    Qmin,A_min=calAQQ2(geo,bed,point)
    Qmean_min=fisher_mean(Qmin)
    N=len(Qmin)
    An=A_min/N
    Dir_remag=[N,point[0],point[1],round(A_min,3),round(An,3)]
    print 'Used direction: ', 'Dec / Inc ', '(',"%.1f" % Dir_remag[1],'/',"%.1f" % Dir_remag[2],')', 'A/n ',"%.3f" % An
else:
    Q1,A1=calAQQ2(geo,bed,point1)
    Q2,A2=calAQQ2(geo,bed,point2)
    Q3,A3=calAQQ2(geo,bed,point3)
    Q4,A4=calAQQ2(geo,bed,point4)
    Qmean_min,Qmin,A_min=minA(geo,bed,Q1,A1)
    Qmean2,Q2,A2=minA(geo,bed,Q2,A2)
    Qmean3,Q3,A3=minA(geo,bed,Q3,A3)
    Qmean4,Q4,A4=minA(geo,bed,Q4,A4)
    if A_min>A2: #looking for the best result coming from the four starting points
        A_min=A2
        Qmean_min=Qmean2
        Qmin=Q2
    if A_min>A3:
        A_min=A3
        Qmean_min=Qmean3
        Qmin=Q3
    if A_min>A4:
        A_min=A4
        Qmean_min=Qmean4
        Qmin=Q4
    N=len(Qmin)
    An=A_min/N
    Dir_remag=[N,round(Qmean_min[0],1),round(Qmean_min[1],1),round(A_min,3),round(An,3)]
    print 'Calculated direction: ', 'Dec / Inc ', '(',"%.1f" % Dir_remag[1],'/',"%.1f" % Dir_remag[2],')', 'A/n ',"%.3f" % An



distance=[] #This is the angle between the BFD and the remagnetization direction, for each site
for dato in Qmin:
    distance_site=ang2point(dato,Qmean_min)
    distance.append(distance_site)

if pregunta_matriz=='y':
    print '\n\nTo be patient! To calculate the matrix spend few min'
    
# In[22]:

paleobed=paleo_dip(tilt,bed,Qmin)


# In[50]:

api=cal_api(geo,bed)

api2=cal_api2(geo,bed)

out_inter=inter(api2,site)

# In[45]:

if pregunta_matriz=='y':    
    mat=cal_matriz(geo,bed)


# In[37]:

#joining the data in a unique list
out_main=[]

for i in range(len(site)):
    site_main=[site[i][0],data[i][1],data[i][2],data[i][3],data[i][4],data[i][5],"%.1f" %tilt[i][0],"%.1f" %tilt[i][1],"%.1f" %Qmin[i][0],"%.1f" %Qmin[i][1],"%.0f" %bed[i][0],"%.2f" %api[i],"%.0f" %paleobed[i][0],"%.0f" %paleobed[i][1],"%.1f" %distance[i]]
    out_main.append(site_main) 
    


# In[25]:

#saving the files
header_main=['Site','Dec_BBC','Inc_BBC','a95', 'DipDir', 'Dip','Dec_ATBC','Inc_ATBC','Dec_BFD','Inc_BFD','Strike','Api_angle','Paleo_strike','Paleo_dip','Remag_dir_angle']
name_main=files[:-4]+'_main'
save_out_file(header_main,out_main,name_main)

header_inter=['Dec','Inc','Site_i','Site_j']
name_inter=files[:-4]+'_inter'
save_out_file(header_inter,out_inter,name_inter)

# In[26]:

if pregunta_remag=='y':
    header_Qmean=['N','Dec','Inc','Asum', 'A/n']
    name_Qmean=files[:-4]+'_Qmean'
    save_out_file(header_Qmean,[Dir_remag],name_Qmean)


# In[36]:

if pregunta_matriz=='y':
    header_matriz=['Dec','Inc','x_cartesian','y_cartesian','x_normalized','y_normalized','A','A/n']
    name_matriz=files[:-4]+'_matrix'
    save_out_file(header_matriz,mat,name_matriz)
    


# In[55]:

time_fin=time()
time_tot=time_fin-time_ini
print '\n\nExecution time: ',round(time_tot,2),'seg'

