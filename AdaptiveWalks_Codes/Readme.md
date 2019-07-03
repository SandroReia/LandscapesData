
1) The code Adaptivewalk-random-hsp90.cpp generates random adaptive walks in the Hsp90 landscape.

As output one has estimates for mean walk length, predictability, mean path divergence, accessibility for each local optimum; and finally fitness values of those local optima.

The data about fitness values and connectivities of the sequences are already embedded in the code, whereas for the estimate
of the mean path divergence the calculation of hamming distance between all pairs of sequences is obtained from the file hamming_distance_tres_colunas.txt

To compile the code:c++ -O3 AdaptiveWalk-random-hsp90.cpp -o AdaptiveWalk-random-hsp90 -lm -lgsl -lgslcblas

To run the code: ./script-random-hsp90

To change the number of adaptive walks just change the script

2) The code Adaptivewalk-prob-hsp90.cpp generates probabilistic adaptive walks in the Hsp90 landscape.

The instructions are the same as the ones for the random version.

3) The code Adaptivewalk-hsp90-randomfreq.cpp generates random adaptive walks in the Hsp90 landscape to calculate the frequency of the mutational pathways produced through the dynamics. Note that the same information is used in 1), but in case one needs a better statistics for the evaluation of the path frequencies, which also warrants that a minimum number of walks is satisfied for every local optimum, the code is more appropriate. The input here is this minimum number of walks terminating at the least visited local optimum.

To compile the code:c++ -O3 AdaptiveWalk-hsp90-randomfreq.cpp -o AdaptiveWalk-hsp90-randomfreq -lm -lgsl -lgslcblas

Tu rum the code: ./script-random-hsp90-freq

4) 



