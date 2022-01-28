# Team-WA2
We are a team from Warsaw University taking part in the Meet EU 2021, a part of the 4EU+ alliance courses. Our topic was "Prediction of TADs".

We are pleased to be sharing with you the outcome of our approach to the problem. We focused on improving the already well-established TopDom algorithm. Since it consists of multiple steps (including signal calculation, smoothing of the signal, peak detection and peak filtering), we reimplemented the original functions in Python and added the possibility for the user to choose their preferred algorithm for each of the algorithm steps. We also chose the set of parameters that to us seem to be giving the most satisfactory results on the data we were provided by the course consultants.

## Download executable files
Download, unpack and run!

## Run our program
Run the PredictTadsWA2.sh script

#### Options
- chromosomes (default: all)
- normalisation (default: normal)
- peak finding function
- smoothing function
- peak filtering function

## Evaluate the results
We have also provided the user with scripts for comparing two outputs to each other or your output to some reference TADs.
#### Overlap score

#### Measure of concordance
