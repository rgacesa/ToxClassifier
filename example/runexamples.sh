#!/bin/bash
echo ' ------------ RUNNING TOX CLASSIFIER EXAMPLE -------------- '
echo '   -> running vectorization; if it crashes, check if paths are fine (see readme)'
echo '../toxClassVectorize.py -I example_input.fa -O _v --verbose 1 --classifynotb 1'
../toxClassVectorize.py -I example_input.fa -O _v --verbose 1 --classifynotb 1
echo '       -> vectorization done!'
echo '   -> running predictor; if it crashes, make sure R and appropriate libraries are installed (see readme)'
echo 'Rscript ../toxClassRClassifier.R --modelsPath ../MLModels/ --vectorsPath ./'
Rscript ../toxClassRClassifier.R --modelsPath ../MLModels/ --vectorsPath ./
echo '       -> DONE, check out.csv, it should look like out_example.csv'

