This program is linked with Calvín et al., (in review). pySCu: a new python code for analyzing remagnetizations directions by means of Small Circle utilities. Computers and geosciences


pySCu has three different modules
 - pySCu_calc.py does the calculations
 - pySCu_draw.py draws different general plots
 - pySCu_draw_labels.py does a plot with the different elements with the name of the site.

The program uses some basic Python libraries as Matplotlib-1.5.3 and Numpy-1.11.2 so will not run on the standard Mac OS and Windows versions of python;
 we recommend either the Anaconda or Canopy installations.
The user is referred to the instructions for PmagPy and Anaconda or Canopy installations in the PmagPy cookbook at: https://earthref.org/PmagPy/. 

Drawing modules are based in PmagPy (Tauxe et al. 2016, G3, DOI: 10.1002/2016GC006307)


Running the programs:


The Python files can be used as executable files but it is recommended to run the program from the command line (against possible errors the command line give us some information about what is happening).
(Command line trick: write the first letters of the program and press 'Tab')


pySCu_calc.py

This is the main programe which does the calculations.
Different workflows are allowed (see workflow_pySC_calc.py) answering two questions (y/n)
Remember that the method allows two symmetric solutions with different (p)ositive or (n)egative inclination (see Symmetry_PmagDir.png)

	What calculations it does?
		- The small circles (SCs)
		- The SCI solution (the remagnetization direction) and its confidence ellipse
		- The paleodips
		- The BFD (best fit direction) and ATBC(after total bedding correction) paleomagnetic directions
		- The SCs intersections
		- The A/n grid

	Input data:
		- spaced delimited text file with header with the next columns (from left to right):
			- Site name
			- Declination of the in situ paleomagnetic direction (declination BBC, before bedding correction)
			- Inclination of the in situ paleomagnetic direction (inclination BBC)
			- Alpha95 of the paleomagnetic direction
			- Kappa (Fisher, 1953) of the paleomagnetic direction
			- Dip direction of the bedding
			- Dip of the bedding
			- Kappa (Fisher, 1953) of the bedding

	Output:
		- A maximum of five files depending of the used workflow. Use the name of the input file appending a last name to differentiate the files:
			- *_main.txt: Main output file with the information relatives to the paleomagnetic directions and the paleobedding
			- *_Ref.txt: The calculated remagnetization direction
			- *_SCIs.txt: The 500 solutions which allows apply the statistic
			- *_matrix.txt: The A/n grid
			- *_inter: the direction of the SCs intersections and the name of the site responsibles of each intersection

pySCu_draw.py

This is the main draw program. It allows us drawing four stereoplots in equal-area projections with:
	- The BBC paleomagnetic directions and the SCs
	- The ATBC paleomagnetic directions and the SCs
	- The best fit direction (BFD) and the SCs
	- The SCs intersections and a contour plot of the A/n value. This plot is optative.

	Input:
		- The program uses the output files from pySCu_calc.py
		- At the starting point, it asks the name of the main file (*_main.txt)
		- Take care: if some of the files are not in the same folder than the *_main.txt file, or were renamed, some errors are possibles in the plot

pySCu_draw_labels.py

It plots an equal-area projection with the SCs, the paleomagnetic directions (BBC, BFD and ATBC) and the name of the site
It can be useful to use in specific points, as to see the relationship between the different sites in a fold, etc.
