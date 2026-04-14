# Metaproteomics_FIT

This is the collection of scripts used for the **"Fecal metaproteomic profiling of leftover FIT samples"** by Avershina et al, 2026

---

## Repository Structure
```sim_input/combined_protein.xslx```: Mock input file with simulated metaproteomic output from the FragPipe

```preprocess.py```: Filtering of the metaproteomics data and formatting input 

```protocol_reproducibility.py```: Finding proteins consistently recovered between stool and FIT samples and calculating their intensities

```samples_variation.py```: PCA analysis of protein intensities; PCoA of Jaccard distances between the samples based on presence/absence of proteins

```check_hydrophobicity.py```: calculate GRAVY scores for the protein sequences

---

## Data Availability

Raw spectra are available from The Federated European Genome-phenome Archive (FEGA): accession xxxxxxxxxxxx

---

## Dependencies

* Python ≥ 3.9
* pandas
* numpy
* scikit-bio
* biopython
* stats
* scanpy
* scipy
* statsmodels
* statannotations
* anndata
* matplotlib
* seaborn

---

## Contact

Ekaterina Avershina<br> 
University of Oslo<br>
Email: ekateria@uio.no
<br>
<br>
Trine B Rounge<br> 
University of Oslo<br>
Email: t.b.rounge@uio.no