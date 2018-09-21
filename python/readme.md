# EEG Analysis

## Data Cleaning
using `mat` format. 
(for scipy.io.loadmat, please use `mat` lower than `v7.3`)

should at least have 3 variables:
- `channels`: m\*n matrix, columns as channels from `ch 1` to `ch m`, rows as raw data points. 
- `times`: n length 1D vector
- `markers`: struct(markername: markervalue)
- `cue_onset`: 1D vector for markers (not used any more)