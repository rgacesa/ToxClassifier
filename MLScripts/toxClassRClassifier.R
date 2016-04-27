# ----------------------------------------------
# ----- TOXIN CLASSIFIER, R script -----------
# ----------------------------------------------
# Does following: 
# 1) loads models 
# 2) loads vectorized input (vectorization is done by python part)
# 3) does prediciton on input
# 4) outputs prediction(s) to file ... or puts it directly to html (hm... ?)
#   - output is in form: sequenceID, prediction STB (Scored Toxbits), prediction TBE (TriBlast Enhanced)
#
#
#

loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

print ('-----------------------------------------------------------------')
print ('-------------- ToxClassifier, predictor script          ---------')
print ('-----------------------------------------------------------------')

# load argparse library (and suppress default output)
suppressPackageStartupMessages(library(argparse))
# ---- here we parse CL arguments --
# create argument parser
parser <- ArgumentParser()
parser$add_argument("--modelsPath",default="/home/ranko/Research/ToxClassifier/MLScripts/test/",help="path to folder with ML models")
parser$add_argument("--vectorsPath", default="/home/ranko/Research/ToxClassifier/testor2/",help="path to folder with input vectors")
parser$add_argument("-O","--output",default="out.csv",help="output file")
parser$add_argument("-L","--log",default="toxClassRClassifier.log",help="log file")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
logline <- " -> parsing CL arguments..."
print(logline)

args <- parser$parse_args()
logline <- "    -> done parsing CL"
write(logline,file=args$log)
print(logline)
logline <- paste("     -> input folder (with models): ",args$modelsPath)
write(logline,file="toxClassR.log",append=TRUE)
print(logline)
logline <- paste("     -> input folder (with vectors): ",args$vectorsPath)
write(logline,file=args$log,append=TRUE)
print(logline)
write(logline,file=args$log,append=TRUE)
print(logline)
logline <- paste("     -> output file: ",args$output)
write(logline,file=args$log,append=TRUE)
print(logline)

# load R caret for ML stuff
logline <- (" -> loading librarires...")
write(logline,file=args$log,append=TRUE)
print(logline)
suppressPackageStartupMessages(library(caret))
logline <- ("    -> libraries loaded")
write(logline,file=args$log,append=TRUE)
print(logline)

# load vector(s) (and check for errors, exit 'gracefully' on error)
print (" -> loading vectorized input and models...")
# set paths and check for errors
vecSTBp <- paste(args$vectorsPath,"_v_stb_vec.csv",sep ="")
vecTBEp <- paste(args$vectorsPath,"_v_tbe_vec.csv",sep = "")
vecBIFp <- paste(args$vectorsPath,"_v_bif_vec.csv",sep = "")

modelSTB_SVMp <- paste(args$modelsPath,"toxScoredSVM40000.RData",sep = "")
modelSTB_GBMp <- paste(args$modelsPath,"toxScoredGBM40000.RData",sep = "")

modelTBEa_SVMp <- paste(args$modelsPath,"ulallSVM53004.RData",sep = "")
modelTBEa_GBMp <- paste(args$modelsPath,"ulallGBM53004.RData",sep = "")

modelTBEb_SVMp <- paste(args$modelsPath,"ul_75_SVM53001.RData",sep = "")
modelBIF_SVMp <- paste(args$modelsPath,"aaBiFreqSVM56540.RData",sep = "")
modelBIF_GBMp <- paste(args$modelsPath,"aaBiFreqGBM56540.RData",sep = "")

modelTBS_SVM <- paste(args$modelsPath,"tbSVM56540.RData",sep="")
modelTBS_GBM <- paste(args$modelsPath,"tbGBM56540.RData",sep="")

if (!file.exists(modelTBS_SVM)) {
    logline <- paste("ERROR: model: file does not exist: ", modelTBS_SVM)
    write(logline,file=args$log,append=TRUE)
    print(logline)    
    stop()
}
if (!file.exists(modelTBS_GBM)) {
    stop(sprintf("ERROR: model: file ( %s ) does not exist", modelTBS_GBM))
}
if (!file.exists(modelSTB_GBMp)) {
    stop(sprintf("ERROR: model: file ( %s ) does not exist", modelSTB_GBMp))
}
if (!file.exists(modelTBEa_SVMp)) {
    stop(sprintf("ERROR: model: file ( %s ) does not exist", modelTBEa_SVMp))
}
if (!file.exists(modelTBEa_GBMp)) {
    stop(sprintf("ERROR: model: file ( %s ) does not exist", modelTBEa_GBMp))
}
if (!file.exists(modelTBEb_SVMp)) {
    stop(sprintf("ERROR: model: file ( %s ) does not exist", modelTBEb_SVMp))
}
if (!file.exists(modelBIF_SVMp)) {
    stop(sprintf("ERROR: model: file ( %s ) does not exist", modelBIF_SVMp))
}
if (!file.exists(modelBIF_GBMp)) {
    stop(sprintf("ERROR: model: file ( %s ) does not exist", modelBIF_GBMp))
}
if (!file.exists(vecSTBp)) {
    stop(sprintf("ERROR: vector: file ( %s ) does not exist", vecSTBp))
}
if (!file.exists(vecTBEp)) {
    stop(sprintf("ERROR: vector: file ( %s ) does not exist", vecTBEp))
}
if (!file.exists(vecBIFp)) {
    stop(sprintf("ERROR: vector: file ( %s ) does not exist", vecBIFp))
}

# no errors, actually load them
#   load modelSTB (requires weird stuff because R sucks)
modelSTB <- loadRData (modelSTB_SVMp)
modelTBE <- loadRData (modelSTB_GBMp)

modelSToxB_SVM <- loadRData(modelSTB_SVMp)
modelSToxB_GBM <- loadRData(modelSTB_GBMp)
modelTBEa_SVM <- loadRData(modelTBEa_SVMp)
modelTBEa_GBM <- loadRData(modelTBEa_GBMp)
modelTBEb_SVM <- loadRData(modelTBEb_SVMp)
modelBIF_SVM <- loadRData(modelBIF_SVMp)
modelBIF_GBM <- loadRData(modelBIF_GBMp)
modelTriOld_SVM <- loadRData(modelTBS_SVM)
modelTriOld_GBM <- loadRData(modelTBS_GBM)

vecSTB <- read.csv(vecSTBp)
vecTBE <- read.csv(vecTBEp)
vecBIF <- read.csv(vecBIFp)
print ("    -> input and models loaded")

# run prediction!
print (" -> running classifiers ...")

if (nrow(vecSTB) < 1) {
    vecSTB <- data.frame(qID = vecBIF$qID)
    predSTB_SVM <- c(NA)
    predSTB_GBM <- c(NA)
} else {
    print ("     -> STB_SVM model (simple tri-blast, support-vector-machine)")    
    suppressMessages(predSTB_SVM <- predict(modelSToxB_SVM,vecSTB))
    print ("     -> STB_GBM model (scored tox-bits, tree boost)")    
    suppressMessages(predSTB_GBM <- predict(modelSToxB_GBM,vecSTB))
}

print ("     -> TBEa_SVM (enhanced tri-blast I, support-vector-machine)")
suppressMessages(predTBEa_SVM <- predict(modelTBEa_SVM,vecTBE))
print ("     -> TBEa_GBM (enhanced tri-blast I, tree boost)")
suppressMessages(predTBEa_GBM <- predict(modelTBEa_GBM,vecTBE))
print ("     -> TBEb_SVM (enhanced tri-blast II, support-vector-machine)")
suppressMessages(predTBEb_SVM <- predict(modelTBEb_SVM,vecTBE))
print ("     -> BIF_SVM (amino-acid bimer frequency, support-vector-machine)")
suppressMessages(predBIF_SVM <- predict(modelBIF_SVM,vecBIF))
print ("     -> BIF_GBM model (amino-acid bimer frequency, tree boost)")
suppressMessages(predBIF_GBM <- predict(modelBIF_GBM,vecBIF))
print ("     -> TBS_SVM (simple tri-blast, support-vector-machine)")
suppressMessages(predTBS_SVM <- predict(modelTriOld_SVM,vecTBE)) # tri blast 'simple'
print ("     -> TBS_GBM (simple tri-blast, tree boost)")
suppressMessages(predTBS_GBM <- predict(modelTriOld_GBM,vecTBE)) # tri blast 'simple'

print (" -> preparing output...")
origOdf <- data.frame(qID=vecBIF$qID)
origOdf$order <- c(1:nrow(origOdf))

predBIF_SVMdf <- data.frame(qID=vecBIF$qID,BIF_SVM=predBIF_SVM)
predBIF_GBMdf <- data.frame(qID=vecBIF$qID,BIF_GBM=predBIF_GBM)
predRes <- merge(predBIF_SVMdf,predBIF_GBMdf,by="qID",sort = FALSE,all.y = TRUE)

resSTB_SVMdf <- data.frame(qID=vecSTB$qID,SToxB_SVM=predSTB_SVM)
predRes <- merge(resSTB_SVMdf,predRes,by="qID",sort = FALSE,all.y = TRUE)

resSTB_GBMdf <- data.frame(qID=vecSTB$qID,SToxB_GBM=predSTB_GBM)
predRes <- merge(resSTB_GBMdf,predRes,by="qID",sort = FALSE,all.y = TRUE)

predTBS_SVMdf <- data.frame(qID=vecTBE$qID,TBSim_SVM=predTBS_SVM)
predRes <- merge(predTBS_SVMdf,predRes,by="qID",sort = FALSE,all.y = TRUE)

predTBS_GBMdf <- data.frame(qID=vecTBE$qID,TBSim_GBM=predTBS_GBM)
predRes <- merge(predTBS_GBMdf,predRes,by="qID",sort = FALSE,all.y = TRUE)

resTBEa_SVMdf <- data.frame(qID=vecTBE$qID,TBEa_SVM=predTBEa_SVM)
predRes <- merge(resTBEa_SVMdf,predRes,by="qID",sort = FALSE,all.y = TRUE)

resTBEa_GBMdf <- data.frame(qID=vecTBE$qID,TBEa_GBM=predTBEa_GBM)
predRes <- merge(resTBEa_GBMdf,predRes,by="qID",sort = FALSE,all.y = TRUE)

resTBEb_SVMdf <- data.frame(qID=vecTBE$qID,TBEb_SVM=predTBEb_SVM)
predRes <- merge(resTBEb_SVMdf,predRes,by="qID",sort = FALSE,all.y = TRUE)

predRes <- merge(origOdf,predRes,by="qID")
predRes <- predRes[order(predRes$order),]
predRes$order <- NULL

write.csv(predRes,args$output)

logline <- ("    ----------------------------------------")
write(logline,file=args$log,append=TRUE)
print(logline)
logline <- ("    -> RUN SUCCESSFUL, ALL DONE !")
write(logline,file=args$log,append=TRUE)
print(logline)
logline <- ("    ----------------------------------------")
write(logline,file=args$log,append=TRUE)
print(logline)
