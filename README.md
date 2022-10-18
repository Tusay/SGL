Replication of results from data:

1. All the data files for the SGL project observations taken on November 6th are stored on the Berkeley archive at:
    http://bldata.berkeley.edu/AGBT21B_999_24/
        Each blc?? folder corresponds to a compute node, each covering a slice of the total frequency range. 
            i.e. the data is split between each blc folder.
        The files can be retrieved using the "download_data.pbs" script, which uses wget. 
        File naming convention is the blc node, guppi (don't ask), the MJD date, target name, target sequence number, then 0000, 0001, or 0002, and finally the suffix .h5
        We named our observations ON_L and OFF_L, for the on-source and off-source pointings, with the letter corresponding to the receiver band (L or S). 
        The three file numbers correspond to frequency resolution. For this project we are looking for narrowband signals, and only care about high resolution files. 
        So, only files labeled ON*0.h5 or OFF**0.h5 are of use. 
        (The normal date listed next to each file in the archive corresponds to the date of transfer to the archive)

2. There are two options for dealing with all the L band and S band data.
    1. Splice together all the nodes into one large file for each 5min observation using splice2.c
        This is what we did for this project. It is no longer recommended. It makes the files huge and difficult to deal with.
        turboSETI works equally well on the unspliced files.
    2. Run turboSETI on each 5min observation in each blc node. 