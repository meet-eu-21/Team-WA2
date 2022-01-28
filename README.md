# Team-WA2
We are a team from Warsaw University taking part in the Meet EU 2021, a part of the 4EU+ alliance courses. Our topic was "Prediction of TADs".

We are pleased to be sharing with you the outcome of our approach to the problem. We focused on improving the already well-established TopDom algorithm. Since it consists of multiple steps (including signal calculation, smoothing of the signal, peak detection and peak filtering), we reimplemented the original functions in Python and added the possibility for the user to choose their preferred algorithm for each of the algorithm steps. We also chose the set of parameters that to us seem to be giving the most satisfactory results on the data we were provided by the course consultants.

## Download executable files
Download the files from *PredictTADsWA2* directory. Unpack and it's ready to run!

## Run our program
Run the *PredictTadsWA2.sh* script with your selected options:

### Options
- ```-c``` chromosomes (default: ```all```)

  chromosomes for which you want the predictions to happen; if more than one, put the list in quotation marks: ```-c "1 2 3"```
  
  if ```-c all```, computations will be performed for 23 chromosomes, ```1 2 3 ... 21 22 X```
  
- ```-n``` normalisation (default: ```normal```)

  an algorithm to be used to normalise the contact matrices 

  ```-n normal``` an algorithm provided by Leopold Carron
  
  ```-n kr``` KR normalisation algorithm; runs longer, but may yield different results

- ```-p``` peak finding function (default: ```find_peaks```, original TopDom uses: ```detect_local_extrema```)

  ```find_min``` ```detect_local_extrema``` ```find_peaks``` ```find_peaks_2```
- ```-s``` smoothing function (default: ```False```, original TopDom uses: ```False```)
  
  ```False``` ```savgol_filter``` ```qspline```
- ```-f``` peak filtering function (default: ```filter_peaks```, original TopDom uses: ```statFilter```)

  ```False``` ```statFilter``` ```filter_peaks```

## Evaluate the results
We have also provided the user with scripts for comparing two outputs to each other or your output to some reference TADs.
### Overlap score
usage: ```Rscript ComputeTopDomOverlapScores.R <1st file> <2nd file> T/F```

The last parameter corresponds to whether you want a full output printed to your console (T) or not (F). The order of the files will only influence the order of your outputs, because the script will always compare the 1st to the 2nd and then the 2nd to the 1st.
### Measure of concordance
usage: ```python3 ComputeMOC.py <1st file> <2nd file>```

The order of the files does not matter.
