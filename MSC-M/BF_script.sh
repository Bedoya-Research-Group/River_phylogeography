#!/bin/sh
 for (( i=1; i<5; i++ )); do
#	python savage_dickey.py 0.01 "M_downmid->up" 1 10 Aguacate/ancestor_up/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_up->downmid" 1 10 Aguacate/ancestor_up/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_mid->down" 1 10 Aguacate/mid_down/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_down->mid" 1 10 Aguacate/mid_down/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_mid->up" 1 10 Aguacate/up_mid/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_up->mid" 1 10 Aguacate/up_mid/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_down->up" 1 10 Cocle/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_up->down" 1 10 Cocle/mcmc"$i".txt
	python savage_dickey.py 0.01 "M_downmid->up" 1 10 Diego/ancestor_up/mcmc"$i".txt
	python savage_dickey.py 0.01 "M_up->downmid" 1 10 Diego/ancestor_up/mcmc"$i".txt
	python savage_dickey.py 0.01 "M_mid->down" 1 10 Diego/mid_down/mcmc"$i".txt
	python savage_dickey.py 0.01 "M_down->mid" 1 10 Diego/mid_down/mcmc"$i".txt
	python savage_dickey.py 0.01 "M_mid->up" 1 10 Diego/up_mid/mcmc"$i".txt
	python savage_dickey.py 0.01 "M_up->mid" 1 10 Diego/up_mid/mcmc"$i".txt	
#	python savage_dickey.py 0.01 "M_downmid->up" 1 10 Paraiso/ancestor_up/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_up->downmid" 1 10 Paraiso/ancestor_up/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_mid->down" 1 10 Paraiso/mid_down/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_down->mid" 1 10 Paraiso/mid_down/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_mid->up" 1 10 Paraiso/up_mid/mcmc"$i".txt
#	python savage_dickey.py 0.01 "M_up->mid" 1 10 Paraiso/up_mid/mcmc"$i".txt	
done
