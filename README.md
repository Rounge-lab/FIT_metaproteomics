# Metaproteomics_FIT

This is the collection of scripts used for the **"Fecal metaproteomics from FIT samples: a feasibility study"** by Avershina et al, 2026

---

## Repository Structure

```input/true_MaxLFQ_nohuman.xslx```: MaxLFQ intensities of FIT and corresponding stool samples from anonimyzed individuals A,B and C and negative controls (only bacterial proteins are reported here)

```input/meta.xslx```: metadata for FIT and corresponding stool samples

```input/simulated_FragPipe_output.xslx```: Mock input file with simulated metaproteomic output from the FragPipe used as input for the ```preprocess.py``` script

```preprocess.py```: Filtering of the metaproteomics data and formatting input 

```check_taxonomy.py```: Finding which bacterial taxa proteins were recovered from

```COG_summary.py```: make a summary of proteins recovered from FIT ans stool grouped by COGs

```samples_variation.py```: PCA analysis of protein intensities; PCoA of Jaccard distances between the samples based on presence/absence of proteins


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