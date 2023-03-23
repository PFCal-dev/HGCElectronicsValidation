# RunMixAndCluster

Will build jets out of signal and pileup nanoaod -like file with CaloParticles and GenParticles

```
a=(`ls *.py`)
for i in ${a[@]}; do
    ../../../../bin/el8_amd64_gcc11/RunMixAndCluster UserCode/HGCElectronicsValidation/bin/${i} &
done
```