# Neoepitope Prediction


Here we will lay out our pipeline for inferring neoepitopes in cancer.

_**Steps**_

* Git clone the repo:
```
git clone https://github.com/tanyaphung/Neoepitope_Prediction.git
```

* Change working directory to Neoepitope_Prediction:
```
cd Neoepitope_Prediction/
```

* Create a conda environment
```
conda env create --name NeoepitopePrediction --file environment.yml
```

* After you create a conda environment, activate it:

```
conda activate NeoepitopePrediction
```

* Download the IEDB tool from http://tools.iedb.org/mhci/download/
 - untar the folder:
 ```
 tar -xzvf IEDB_MHC_I-2.19.1.tar.gz
 ```
 - configure the tool and exit out of the directory 
 ```
 cd mhc_i/
 ./configure
 cd ..
 ```
 
 - You will get a message like this after configuring the tool:
 ```
 All prerequisites found!
Copying the standalone-specific netMHCcons template into place
IEDB MHC class I binding prediction tools successfully installed!
Use the command 'python src/predict_binding.py' to get started
 ```

* Add peptides to respective folder. 
 - Note that the file format has to be `"TCGA-" + patient +"_Varscan_variants_filter.pass."+ str(num) +".peptide"`.
* Add hla to the hla folder. Then delete the file `all_hlas_here.txt`. 
```
rm -f hla/all_hlas_here.txt
```

* Submit job on cluster:
 - Modify the script `submit.sh` to put in your information
 - Note: because of some issue in memory, right now this will work on ASU agave cluster. 
