# Python Analysis Code

Python scripts and Jupyter notebooks used to analyze simulation results. 

* `histograms.py` -- Simulations recorded raw avalanches. This python script takes the raw list of data and çreates a new file containing the log-binned histograms. 
* `meanCluster.py` -- Takes in the log-binned histograms. Returns the mean cluster size (for a given parameter set). 
* `correlationLength.py` -- Takes in the log-binned histograms. Returns the correlation length (for a given parameter set). 
* `avg_sizevsDur.py` -- Takes in raw list of avalanche sizes and durations  Returns the average size of avalanches for a given duration. 
* `mm_simul_distr.py` -- Simulates the Manna model with the multinomial distribution algorithm. Written in Python for small, quick simulations. 