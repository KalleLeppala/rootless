# rootless

This repository contains the code and some simulated data used in the article “Notes on *f<sub>4</sub>*-ratio estimation”, Kalle Leppälä.

The simulation of a bidirectional gene flow scenario and testing the bidirectional admixture proportion estimators are described in the files bidirectional.slim, bidirectional.py and bidirectional.R. The files bidirectional.py and bidirectional.R also contain (commented out) PLINK commands: they are meant to be read and executed manually. The folder merged contains the precomputed blockwise *f<sub>2</sub>*-statistics constructed early in bidirectional.R, the full data take a while to create and are too big for storing here.

The benchmarking of the enhanced statistic *f<sub>B</sub>* using a simulated cichlid data set is described in cichlids.R.

The bear admixture example is described in bears.R.