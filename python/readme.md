# EEG Analysis

## Data Cleaning
using `mat` format. 
(for scipy.io.loadmat, please use `mat` lower than `v7.3`)

should at least have 3 variables:
- `channels`: m\*n matrix, columns as channels from `ch 1` to `ch m`, rows as raw data points. 
- `times`: n length 1D vector
- `markers`: struct(markername: markervalue)
- `cue_onset`: 1D vector for markers (not used any more)

## TODO

### Entrain
- [x] average frequency-time domain output
    - [x] grating (tf_domain)
    - [x] entrain (tf_domain)
    - [x] preview (preview/tf_domain)
- [ ] latency time
    - [x] grating (latency)
    - [x] entrian (latency)
    - [ ] preview (preview/bandpower)
- [ ] auc in grating periods and inter-stimuli intervals (auc)
- [ ] time points of all events
    - [ ] threshold (event_times)
    - [ ] range (event_times)
- [ ] raw data
    - [ ] raw data preview (preview/raw)
    - [ ] bandpass filted preview (preview/bandpass)