#FastQC
export PATH=/opt/exp_soft/qiime/packages/other/FastQC/:$PATH
# QIIME specific
export PYTHONPATH=/opt/exp_soft/qiime/lib/python2.7/site-packages/:$PYTHONPATH
## usearch would also be linked from her
export PATH=/opt/exp_soft/qiime/bin/:$PATH
export PATH=/opt/exp_soft/qiime/python-2.7.3/bin:$PATH
export RDP_JAR_PATH=/opt/exp_soft/qiime/app/rdp_classifier_2.2/rdp_classifier-2.2.jar
export PATH=/opt/exp_soft/qiime/R-3.0.2/bin:$PATH
export PATH=$PATH:/opt/exp_soft/qiime/app/ChimeraSlayer
export PYRO_LOOKUP_FILE=/opt/exp_soft/qiime/packages/extras/AmpliconNoiseV1.27/Data/LookUp_E123.dat
export SEQ_LOOKUP_FILE=/opt/exp_soft/qiime/packages/extras/AmpliconNoiseV1.27/Data/Tran.dat
export SOURCETRACKER_PATH=/opt/exp_soft/qiime/packages/extras/sourcetracker-0.9.5
export BLASTMAT=/opt/exp_soft/qiime/packages/extras/blast-2.2.22/data
export DONT_USE_MPI=1

# Additional software
export PATH=$PATH:/opt/exp_soft/qiime/packages/other
export PATH=/opt/exp_soft/qiime/packages/other/ImageMagick-6.9.3-5/install/bin/:$PATH
export PATH=$PATH:/opt/exp_soft/seqtk

# Tutorial specific code
export PATH=$PATH:/scratch/DB/bio/training/16SrRNA/16SrRNA-hex-tutorial/src
export PATH=$PATH:/scratch/DB/bio/training/16SrRNA/16SrRNA-hex-tutorial/src/fastqc_combine
