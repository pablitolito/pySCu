# -*- coding: utf-8 -*-

import csv
import numpy as np
import sys
import pmag as pmag
import pyscu_libs as scu
rad=np.pi/180
deg=180/np.pi
from time import time

def main():

	#The input file
	files = raw_input("\nFile name: ")
	try:
		file = open(files)
	except:
		print "Try again (do not forget the extension of the file)"
		exit()


	data,site,geo,bed_d,bed_s,bed_pole,N_sites=scu.saveInputFile(files)   #Saving in diferent list the data
	tilt=scu.tilt_rot(geo,bed_s)                   #Calculating TILT directions

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

	time_ini=time()

	#Calculating BFD and A/N for a manual input remagnetization direction.
	if pregunta_remag=='n':
		ref=[DecInp, IncInp]
		BFD,A_min=scu.calAQQ2(geo,bed_s,ref)
		Qmean_min=scu.fisher_mean(BFD)
		N=len(BFD)
		An=A_min/N
		Dir_remag=[N,ref[0],ref[1],round(A_min,3),round(An,3)]
		print 'Used direction: ', 'Dec / Inc ', '(',"%.1f" % Dir_remag[1],'/',"%.1f" % Dir_remag[2],')', 'A/n ',"%.3f" % An
	else: #Calculating 200 SCI_solutions, their mean (the remagnetization direction) and the rest of parameters
		SCIs=[]
		pgeo=scu.para_dir(geo)
		pbed_pole=scu.para_dir(bed_pole)
		aaa=500 #Change this value for a different value of SCI_solutions calculation
		w,cont=0,aaa/10
		print '\nplease, be patient... calculating',aaa,'SCI_solutions'
		for l in range (aaa): 
			w+=1
			dec=np.random.randint(0,359,1)
			inc=np.random.randint(1,89,1)
			if pregunta_inc=='n': inc=inc*(-1)
			point=[dec[0],inc[0]]
			pgeo_u=scu.selec_para_geo(pgeo)
			pbed_u=scu.selec_para_pole(pbed_pole)
			pQ=scu.p_calAQQ2(pgeo_u, pbed_u, point)
			#print pgeo_u
			#print pbed_u
			#print pQ
			Qmean_min=scu.p_minA(pgeo_u, pbed_u,pQ)
			SCIs.append(Qmean_min)
			#print SCIs
			if w==cont:
				print cont
				cont+=aaa/10
			
		Qend=pmag.dokent(SCIs,1.)
		#print Qend
		ref=[Qend['dec'], Qend['inc']]
		#print ref
		BFD,A_min=[],[]
		BFD,A_min=scu.calAQQ2(geo,bed_s,ref)
		An=A_min/N_sites
		




	distance=[] #This is the angle between the BFD and the remagnetization direction, for each site
	for dato in BFD:
		distance_site=scu.ang2point(dato,ref)
		distance.append(distance_site)


	paleobed=scu.paleo_dip(tilt,bed_s,BFD)


	api=scu.cal_api(geo,bed_s)
	api2=scu.cal_api2(geo,bed_s)

	out_inter=scu.inter(api2,site)
	
	mr=len(out_inter)
	mp=N_sites*(N_sites-1.)/2.
	mr_mp=round(mr/mp,2)
	
	

	#joining the data in a unique list
	out_main=[]
		
	Dir_remag=[N_sites,round(Qend['dec'],1),round(Qend['inc'],1),round(Qend['Eta'],1),round(Qend['Zeta'],1),
	round(Qend['Edec'],1),round(Qend['Einc'],1),round(Qend['Zdec'],1),round(Qend['Zinc'],1),round(An,3),round(A_min,3),mr_mp,aaa]
	print '\nKent mean remagnetization direction (Kent, 1982; Tauxe et al., 1991)', '\nDec / Inc: ', "%.1f" % Dir_remag[1],"/","%.1f" % Dir_remag[2]
	print 'A/n: ',"%.3f" % An,'mr/mp: ',"%.2f" % mr_mp
	print 'Eta_95, dec, inc:', "%.1f" % Dir_remag[3], "%.1f" % Dir_remag[5], "%.1f" % Dir_remag[6]
	print 'Zeta_95, dec, inc', "%.1f" % Dir_remag[4], "%.1f" % Dir_remag[7], "%.1f" % Dir_remag[8] 	
	
	for i in range(len(site)):
		site_main=[site[i][0],data[i][1],data[i][2],data[i][3],data[i][4],data[i][5],data[i][6],data[i][7],"%.1f" %tilt[i][0],"%.1f" %tilt[i][1],"%.1f" %BFD[i][0],"%.1f" %BFD[i][1],"%.0f" %bed_s[i][0],"%.2f" %api[i],"%.0f" %((paleobed[i][0]+90)%360),"%.0f" %paleobed[i][1],"%.1f" %distance[i]]
		out_main.append(site_main) 
		

	#saving the files
	header_main=['Site','Dec_BBC','Inc_BBC','a95', 'k','DipDir', 'Dip', 'k_bed','Dec_ATBC','Inc_ATBC','Dec_BFD','Inc_BFD','SC_Strike','SC_Api_angle','Paleo_DipDir','Paleo_dip','Ref_BFD_angle']
	name_main=files[:-4]+'_main'
	scu.save_out_file(header_main,out_main,name_main)

	header_inter=['Dec','Inc','Site_i','Site_j']
	name_inter=files[:-4]+'_inter'
	scu.save_out_file(header_inter,out_inter,name_inter)

	# In[26]:

	if pregunta_remag=='y':
		header_Ref=['N','Dec','Inc','Eta','Zeta','Dec_E','Inc_E','Dec_Z','Inc_Z','A/n','Asum','mr_mp','number of SCI_solutions']
		name_Ref=files[:-4]+'_Ref'
		scu.save_out_file(header_Ref,[Dir_remag],name_Ref)
		header_SCIs=['Dec','Inc']
		name_SCIs=files[:-4]+'_SCIs'
		scu.save_out_file(header_SCIs,SCIs,name_SCIs)

	# In[36]:

	if pregunta_matriz=='y':
		print '\n\nTo be patient! To calculate the matrix spend few min'
		mat=scu.cal_matriz(geo,bed_s)
		header_matriz=['Dec','Inc','x_eqarea','y_eqarea','x_eq_normalized','y_eq_normalized','A','A/n']
		name_matriz=files[:-4]+'_matrix'
		scu.save_out_file(header_matriz,mat,name_matriz)
		


	# In[55]:

	time_fin=time()
	time_tot=time_fin-time_ini
	print '\n\nExecution time: ',round(time_tot,2),'seg'
	

main()

