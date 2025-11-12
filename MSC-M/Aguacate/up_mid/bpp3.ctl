		seed = -1
		
	seqfile =  ./MSCM_Aguacate.phy
	Imapfile = ./imap.txt
	outfile =  ./out3.txt
	mcmcfile = ./mcmc3.txt

speciesdelimitation = 0
		speciestree = 0

	 species&tree = 3		up down mid
							 1 5 3
							((down,mid), up);

		phase = 1 1 1
		usedata = 1
		nloci = 199
		  
	  cleandata = 0

      migprior = 1 10
      migration = 2 
					up mid
					mid up
							  
	 thetaprior = 3 0.0004	e 	
	   tauprior = 3 0.0008  
	   
       finetune = 1: 0.01	0.02	0.03	0.04	0.05	0.01	0.01
       
          print = 1	0	0	0	0
         burnin = 20000
       sampfreq = 2
        nsample = 200000
        threads = 3
