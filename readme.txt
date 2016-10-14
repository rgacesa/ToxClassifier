ToxClassifier, developed by:
 - Ranko Gacesa, King's College London under supervision of Dr. Paul Long, King's College London
 - last update: 28/04/2016

------------------------------------------------------------------------------------
ToxClassifier code
------------------------------------------------------------------------------------
This is 'raw' code for toxic protein identification and classification by machine learning methods. 
See http://bioserv7.bioinfo.pbf.hr/ToxClassifier/ and appropriate article (listed on web page) for details

------------------------------------------------------------------------------------
To use offline classifier, do following steps:
------------------------------------------------------------------------------------
0) keep following in mind: 
- use online variant (http://bioserv7.bioinfo.pbf.hr/ToxClassifier/) if possible; it includes extra code for handling problems with input sequences and for 'extra intelligence' to handle problematic sequences which might get dubious classificaiton by offline variant; it also provides user friendly interface, output and help with manual annotation

- code was developed under Ubuntu Linux 12.04, it was tested under linux and was never intended to work or tested on different platforms; it might work on different distributions of Windows, Mac or even Windows, but was never tested on those
    - it is however based on Python and R, both of which are reasonably tolerant towards different platforms, so in principle it should work on different linux distributions and might even work on Mac or Windows with certain adjustments
    - also note some of classifiers can get memory hungry, especially if run with lots of (10s or 100s of thousands of sequences); code was developed on 64 Gb memory machine and runs smoothly on it, it was not tested for memory requirements and while it should work with much less, there is no guarantee it will
- this is open source code distributed under GNU GENERAL PUBLIC LICENSE, it is distributed 'as is' without any guarantee it will actually work on anything except for test system where it was developed; user accepts responsability for any damage or problems it might cause when run!
- we can at best provide very limited support for running this code locally; contact ranko.gacesa@kcl.ac.uk with questions and issues; for running it on large datasets it might be easier for us to run the sample then provide support for installing it on different system
- please cite the appropriate article (see http://bioserv7.bioinfo.pbf.hr/ToxClassifier/) if you used this code or its web-service variant

Actually running the code:
--------------------------

1) prepare sequences
 - make sure all input is in FASTA format, ideally with no suspicious, duplicate or otherwise problematic headers
  - safest way here is to just name them by ID or sequentially (example: >gi_xxx); if headers are problematic, classifier might misbehave as it doesn't have very sophisticated way to deal with this and vectorization vs prediction scripts might get in conflict with header names (web variant has extra code to handle this, this one DOES NOT)
  - put all sequences into single file, it is also good idea to put it into new folder as code generates some random temporary files and DOES NOT CLEAN IT UP (again web variant does)
  - scripts are configured to work in one folder inside ToxClassifierGit folder (for example ToxClassifierGit/example or ToxClassifierGit/analysis1... if not run from there, paths WILL NEED TO BE ADJUSTED!
  - technically there is no limit to how many sequences can be in input file; most of running time is taken by loading models rather then classification, so it is much faster to classify sequences in bulk then one by one; however very large datasets can eat up a lot of memory and it will crash if whole thing cannot fit into memory; also vectorization of very large datasets can take a while

2) run vectorization script (toxClassVectorize.py); 
  --help option will display help
  - make sure folders and paths are set properly (see Warnings/notes and Example); also make sure permissions are set correctly (hmmsearch and blastp need to be executable, 
toxClassVectorize.py needs to be executable or run as "python toxClassVectorize.py")
  - use --classifynotb 1 -O _v 
     -> this makes sure vectorization of sequences with no toxbits is done properly (--classifynotb 0 is for making training sets!) and
output is named appropriately for input into classifier (vectors should be named: _v_bif_vec.csv, _v_stb_vec.csv, _v_tbe_vec.csv)

3) from the the same folder, use Rscript to run toxClassRClassifier.R; ideally this is one folder inside, so command should be 
Rscript ../toxClassRClassifier.R

4) output is saved in same folder, as out.csv
 - it is csv file which lists result of each classifier for each input sequence
 - in general, higher total score means sequence is more likely to be toxin; however, if there are no close toxins at all, it should be treated as suspicous at best and if there are very close toxins with score higher then close non-toxins, it is likely to actually be toxin, no matter the classification score; these two checks are not automatic in offline variant (__tb__tmp_neg1.bres, __tb__tmp_neg2.bres, __tb__tmp_pos.bres) are output of this part of analysis and codes for parsing and making sense of it are not included 
 -see http://bioserv7.bioinfo.pbf.hr/ToxClassifier/ and appropriate article (listed on web page) for more detail and interpretation of results

--------------------------------------------------------------------------------------
To run example (it should run succesfully if everything is configured properly and 
all R libraries are in place)
--------------------------------------------------------------------------------------
1) go to example folder
2) run:
../toxClassVectorize.py -I example_input.fa -O _v --verbose 1 --classifynotb 1
3) run:
Rscript ../toxClassRClassifier.R --modelsPath ../MLModels/ --vectorsPath ./
4) output should appear in ./out.csv, it should look something like: 

OR 
1) go to example folder
2) run runexamples.sh (make sure it is executable)

"","qID","TBEb_SVM","TBEa_GBM","TBEa_SVM","TBSim_GBM","TBSim_SVM","SToxB_GBM","SToxB_SVM","BIF_SVM","BIF_GBM"
"1",">seq1","1","1","1","1","1","1","1","1","1"
"2",">seq2","1","1","1","1","1","1","1","1","1"
"3",">seq3","0","0","0","0","0","0","0","0","0"
"4",">seq4","0","0","0","0","0","0","0","0","0"

-> output is csv format with classifier results (0 = not toxin, 1 = toxin) for each input sequence

Warnings/notes:

A) Vectorization script (toxClassVectorize.py)
- this is python script, written for python2.7, UBUNTU linux 12.04; it was not tested under other environments and might (or not) work
- make sure all paths are in order (default ones WILL NOT WORK! and have to be entered as command line options
or as changed in code - edit following file: 
  - if something is wrong with paths, code MIGHT NOT CRASH, but produce junk instead! Check output from classifier, 
BLAST and HMMER will indicate if they cannot find database(s). 
  - scripts are configured to run from this folder (example) or another folder one level inside ToxClassifier root, paths will need adjusting if run from different location
  - tools folder includes hmmsearch and phmmer; dbs folder includes databases for BLAST and HMMER used for vectorization; these are all absolutely required to do vectorization; 

  to adjust paths, edit toxClassVectorize.py  
      - line 50 - 59: replace paths with appropriate ones on the local machine

B) prediction script (toxClassRClassifier.R)
  - this is R script, it requires R to run, it also requires R caret package 
and selection of other packages connected to R caret; make sure they are all installed
and working correctly   
  - also check paths if models cannot be found
  - machine learning models are in MLModels folder, classifier is utterly useless without them and will not work
