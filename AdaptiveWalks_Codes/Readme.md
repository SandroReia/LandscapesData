
The code Adaptivewalk-random-hsp90.cpp generates random adaptive walks in the Hsp90 landscape.

As output one has estimates for mean walk length, predictability, mean path divergence, accessibility for each local optimum; and finally fitness values of those local optima.

The data about fitness values and connectivities of the sequences are already embedded in the code, whereas for the estimate
of the mean path divergence the calculation of hamming distance between all pairs of sequences is obtained from the file hamming_distance_tres_colunas.txt

To compile the code:c++ -O3 AdaptiveWalk-random-hsp90.cpp -o AdaptiveWalk-random-hsp90 -lm -lgsl -lgslcblas

To run the code: ./script-random-hsp90

To change the number of adaptive walks just change the script


