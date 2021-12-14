#!/bin/bash

for i in {1..5};
do
    echo $i
    python2 to_dcd.py '/salilab/park2/dibyendu/veena_transfer/rigid/very_good_models/all_models_2/very_good_models/ZZ_cluster'$i'_/' './cluster'$i'/cluster.'$i'.sample_A.txt' './cluster'$i'/cluster'$i'_A.dcd' 30
    python2 to_dcd.py '/salilab/park2/dibyendu/veena_transfer/rigid/very_good_models/all_models_2/very_good_models/ZZ_cluster'$i'_/' './cluster'$i'/cluster.'$i'.sample_B.txt' './cluster'$i'/cluster'$i'_B.dcd' 30

done
