"""Module for an automatic download of DoRothEA interactions"""

import os
import datetime

def extraction(organism: str="human", confidence: str="AB", directory_output: str="./"):
    """ 
    Store in a SIF file the subnetwork from DoRothEA database corresponding to interactions having 
    - the choosen confidence level, 
    - about a choosen organism.
    
    INPUT
     * organism: string that can be human or mouse
     * confidence: string in the set A, AB (default), ABC, ABCD, ABCDE
     * output directory: the current one by default
        
    OUTPUT
     * SIF file in the directory_output, 
     under the format "dorothea_<confidence>_<organism>_YYYYMMDD.sif" 
     - with <confidence> the confidence levels given in argument
     - with <organism> the organism given in argument
    """
    
    assert organism == 'human' or organism == 'mouse', f"organism must be human or mouse"
    
    date = datetime.datetime.now()
    
    import rpy2.robjects as robjects
    import rpy2.robjects.packages as rpackages

    if confidence == "A":
        confidenceR = '"A"'
    elif confidence == "AB":
        confidenceR = '"A","B"'
    elif confidence == "ABC":
        confidenceR = '"A","B","C"'
    elif confidence == "ABCD":
        confidenceR = '"A","B","C","D"'
    elif confidence == "ABCDE":
        confidenceR = '"A","B","C","D","E"'
    else:
        raise InputError("Incorrect argument: confidence for edges can be A, AB, ABC, ABCD, ABCDE")

    #robjects.r('''
    #    if (!requireNamespace("BiocManager", quietly = TRUE))
    #        install.packages("BiocManager")
    #    BiocManager::install("dorothea")
    #''')
    
    robjects.r('''
        library(dorothea)
        subset_dth = dorothea_{0}[dorothea_{0}$confidence %in% c({1}), ]
        '''.format('hs' if organism=='human' else 'mm', confidenceR))
    robjects.r('''
        df = data.frame(source = subset_dth$tf,
                        sign = subset_dth$mor,
                        target = subset_dth$target)
        write.table(df, file="{0}dorothea_{2}_{3}_{1}.sif", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
        '''.format(directory_output, date.strftime("%Y%m%d"), confidence, organism))