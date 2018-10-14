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

#### 2018-09-30 (discussion)
- [ ] **artifacts remove**
- [ ] entrain quantification
    - [ ] grating power in 2 trials
- [ ] band power curve of 5-15hz and gamma band
- [ ] marker bias, recalculate the auc for clustering
- [ ] relative latency
- [ ] brain areas of electrodes


#### 2018-09-29
- [x] raw waveform averaged across days
    - [x] split data
    - [x] latency -> offset
    - [x] ERP plot
    - [x] normalized total power
- [x] clustering analysis
    - [x] get auc from (-1, 1.5)
    - [x] kmeans
    - [x] DBSCAN
- [ ] phase analysis, biological meaning
    - [ ] phase lock
    - [x] ITPC
    - [ ] wITPCz
    - [ ] PAC

#### 2018-09-26 (discussion)
- [x] Chen Zhou tfdomain plots
    - [x] 10 sec adaptation
    - [x] plots with position name
- [ ] channel classification into 4 categories
    - [x] AUC
    - [x] PCA
    - [x] k-means

- [x] bandpower quantification (by AUC)

- [x] mice date
- [x] phase analysis

#### 2018-09-21 (discussion)
- [x] average frequency-time domain output
    - [x] grating (tf_domain)
    - [x] entrain (tf_domain)
    - [x] preview (preview/tf_domain)
- [x] latency time
    - [x] grating (latency)
    - [x] entrian (latency)
    - [x] preview (preview/bandpower)
- [x] auc in grating periods and inter-stimuli intervals (auc)
- [ ] time points of all events
    - [x] threshold (event_times)
    - [ ] range (event_times)
- [x] raw data
    - [x] raw data preview (preview/raw)
    - [x] bandpass filted preview (preview/bandpass)
- [ ] entrain average
