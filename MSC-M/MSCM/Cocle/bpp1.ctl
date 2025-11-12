		seed = -1
		
	seqfile =  ./MSCM_Cocle.phy
	Imapfile = ./imap.txt
	outfile =  ./out1.txt
	mcmcfile = ./mcmc1.txt

speciesdelimitation = 0
		speciestree = 0

	 species&tree = 2		up down
							3 6
							(up, down);

		phase = 1 1
		usedata = 1
		nloci = 199
		  
	  cleandata = 0

      migprior = 1 10
      migration = 2 
					up down
					down up
							  
	  
	  	 thetaprior = 3 0.0004	e 	
	   tauprior = 3 0.0008  
	   
       finetune = 1: 0.01	0.02	0.03	0.04	0.05	0.01	0.01
       
          print = 1	0	0	0	0
         burnin = 20000
       sampfreq = 2
        nsample = 200000
        threads = 4
