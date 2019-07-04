This folder contains the codes that generate the Gb1 input files for the Adaptative Walks Codes. 

In a few words, the two codes presented here gets the data from the Gb1 empirical landscape and split it into files that are used as input in the Adaptative Walks Codes. The Gb1 empirical landscape files, **elife-16965-supp1-v4.xlsx** and **elife-16965-supp2-v4.xlsx**, were obtained from N. C. Wu, L. Dai, C. A. Olson, J. O. Lloyd-Smith, and R. Sun, Elife 5, e16965 (2016).

To generate the input files, all the files presented here have to be downloaded to the same folder and:

1 - run the python script **preparing_Elife_1.py**. One can do that by typing in the terminal (or powershell): 
> python preparing_Elife_1.py

2 - run the fortran code **preparing_Elife_2.f90**. This can be done by typing in the terminal (or powershell): 
> gfortran preparing_Elife_2.f90 -o elife

and then
> ./elife

As a result of the actions (1) and (2) described above, four files will be created. These files are used as input files for the Adaptative Walks Codes.
