
1) The code Adaptivewalk-random-hsp90.cpp generates random adaptive walks in the Hsp90 landscape. The empirical data from the Hsp90 fitness landscape can be obtained from C. Bank et at. PNAS 113, 14085 (2016).

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

Tu run the code: ./script-random-hsp90-freq

4) The code Adaptivewalk-hsp90-probfreq.cpp generates probabilistic adaptive walks in the Hsp90 landscape to calculate the frequency of the mutational pathways produced through the dynamics. So, the remaining information is exactly the same as in 3). 

To compile the code:c++ -O3 AdaptiveWalk-hsp90-probfreq.cpp -o AdaptiveWalk-hsp90-probfreq -lm -lgsl -lgslcblas

Tu run the code: ./script-probabilistic-hsp90-freq

5) The code Adaptivewalk-GB1-randomfreq.cpp generates random adaptive walks in the GB1 fitness landscape to calculate the frequency of the mutational pathways produced through the dynamics, but also mean walk length, predictability, mean path divergence and accessibility for each local optimum. Likewise, the code warrants that a minimum number of walks is satisfied for every local optimum, which is provided in the script. This is bit tricky, as one of the local optimum of the GB1 landscape is poorly visited through the walks starting at the wild type sequence. All the information about the GB1 landscape is provided by the processed information and contained in the files elife_seq_number.txt, elife_sequence_degree.txt, elife_sequence_fitness.txt and elife_sequence_neighbors_correta.txt. The latter one could not be uploaded (430 Mb) and we can provide upon request. Anyhow, all the information (not processed to be used by the code) about the GB1 landscape is given in Wu et al. Elife 5, e16965 (2016).

To compile the code:c++ -O3 AdaptiveWalk-GB1-randomfreq.cpp -o AdaptiveWalk-GB1-randomfreq -lm -lgsl -lgslcblas

Tu run the code: ./script-GB1-randomfreq





