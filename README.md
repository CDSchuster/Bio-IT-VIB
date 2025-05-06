# Bio-IT-VIB
------------


Code, results and report for the Bio-IT task I was assigned for my second interview for the postdoc position at the De Preter lab in VIB-Ghent.

## Folder structure

1. `Data`: contains the 2 csv files that were provided to perform deconvolution analysis
2. `Outputs`: csv files produced at different stages of the analysis carried out in this project
    * `bulk_samples_deconv` files: the output of meth_atlas analysis
    * `sample_clusters.csv`: the clusters determined in `Data_exploration.ipynb` using PCA and Pearson correlation
    * `cpg_clus_colors_csv`: the colors used to cluster the CpGs in the heatmaps
    * `treatXX_annot.csv` files: the differentially methylated CpGs and the associated genes that were obtained using *R-limma*
    * `Figures`: a folder with all the plots produced in the jupyter notebooks and that were included in the report
3. `Scripts`: All scripts used throughout this project

## Scripts order

Bear in mind, that in order to run the scripts, you will need to change the paths inside the scripts whenever data is loaded or saved

1. `install_conda_env.sh`: this will install the conda environment necessary to run *meth_atlas* from the file `env.yaml`
2. `run_meth_analysis.sh`: runs deconvolution analysis with *meth_atlas* on the 2 files in `Data` folder
3. `Data_exploration.ipynb`: a python notebook where correlation between samples and between CpGs is explored in search of possible clusters
4. `Diff_met_analysis.R`: an R script that performs differential methylation analysis using *R-limma* and returns the DMGs (differentially methylated CpGs) along with annotated genes
5. `Post_DMA.ipynb`: a python notebook to see correlation between the CpG clusters that were determined previously in `Data_exploration.ipynb` with the DMGs. Afterwards, enrichment analysis with *gseapy* is performed on the DMGs to have a better understanding of the biology behind