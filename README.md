# pySCu
Paleomagnetic tool to apply Small Circle methods

This program is linked to Calvin et al. (2017). pySCu: A new Python code for analyzing remagnetization directions using small circle utilities. Computers & Geosciences, 109(March), 32–42. https://doi.org/10.1016/j.cageo.2017.07.002.

### pySCu consists of three modules:

    pySCu_calc.py: Performs the calculations.
    pySCu_draw.py: Generates various general plots.
    pySCu_draw_labels.py: Creates a plot with different elements labeled by site name.


A previous installation of Pmagpy (https://pmagpy.github.io/PmagPy-docs/intro.html) will ensure that all the libraries needed by the programme are installed.

The drawing modules are based on PmagPy (Tauxe et al., 2016, Geochemistry, Geophysics, Geosystems, DOI: 10.1002/2016GC006307).


### Running the programs:

### pySCu_calc.py

This is the main program that performs the calculations. It supports various workflows (see workflow_pySC_calc.py) by answering two questions (y/n). Note that the method allows for two symmetric solutions with either positive (p) or negative (n) inclination (see Symmetry_PmagDir.png).

	What calculations it does?
		- The small circles (SCs)
		- The SCI solution (the remagnetization direction) and its confidence ellipse
		- The paleodips
		- The BFD (best fit direction) and ATBC(after total bedding correction) paleomagnetic directions
		- The SCs intersections
		- The A/n grid

	Input data:
		- A spaced delimited text file with a header contaning the following columns (from left to right):
			- Site name
			- Declination of the in situ paleomagnetic direction (BBC, before bedding correction)
			- Inclination of the in situ paleomagnetic direction (BBC)
			- Alpha95 of the paleomagnetic direction
			- Kappa (Fisher, 1953) of the paleomagnetic direction
			- Dip direction of the bedding
			- Dip of the bedding
			- Kappa (Fisher, 1953) of the bedding

	Output:
		- Up to five files depending of workflow used. These files are named after the input file, with an additional suffix to differentiate them:
			- *_main.txt: Main output file with the information relatives to the paleomagnetic directions and the paleobedding
			- *_Ref.txt: The calculated remagnetization direction
			- *_SCIs.txt: The 500 solutions which allows apply the statistic
			- *_matrix.txt: The A/n grid
			- *_inter: the direction of the SCs intersections and the name of the site responsibles of each intersection

### pySCu_draw.py

This is the primary plotting program. It creates four stereographic projections in equal-area plots showing:
	- The BBC paleomagnetic directions and the SCs
	- The ATBC paleomagnetic directions and the SCs
	- The best fit direction (BFD) and the SCs
	- The SCs intersections and a contour plot of the A/n value. This plot is optative.

	Input:
		- The program uses the output files from pySCu_calc.py
		- At the starting point, it asks the name of the main file (*_main.txt)
		- Take care: if some of the files are not in the same folder than the *_main.txt file, or were renamed, some errors are possibles in the plot

### pySCu_draw_labels.py

This module creates an equal-area projection plot displaying the SCs, paleomagnetic directions (BBC, BFD, and ATBC), and site names. It is useful for specific cases, such as examining relationships between different sites in a fold, etc.
