#!/bin/sh
for (( i=1; i<5; i++ )); do
	python ./savage_dickey.py 0.01 "M_Diego->Aguacate" 1 10 ./run$i/mcmc1.txt
	python ./savage_dickey.py 0.01 "M_Aguacate->Diego" 1 10 ./run$i/mcmc1.txt
	python ./savage_dickey.py 0.01 "M_Diego->Cocle" 1 10 ./run$i/mcmc2.txt
	python ./savage_dickey.py 0.01 "M_Cocle->Diego" 1 10 ./run$i/mcmc2.txt
	python ./savage_dickey.py 0.01 "M_Diego->Piedras" 1 10 ./run$i/mcmc3.txt
	python ./savage_dickey.py 0.01 "M_Piedras->Diego" 1 10 ./run$i/mcmc3.txt
	python ./savage_dickey.py 0.01 "M_Aguacate->Piedras" 1 10 ./run$i/mcmc4.txt
	python ./savage_dickey.py 0.01 "M_Piedras->Aguacate" 1 10 ./run$i/mcmc4.txt
	python ./savage_dickey.py 0.01 "M_Aguacate->Cocle" 1 10 ./run$i/mcmc5.txt
	python ./savage_dickey.py 0.01 "M_Cocle->Aguacate" 1 10 ./run$i/mcmc5.txt
	python ./savage_dickey.py 0.01 "M_Cocle->Piedras" 1 10 ./run$i/mcmc6.txt
	python ./savage_dickey.py 0.01 "M_Piedras->Cocle" 1 10 ./run$i/mcmc6.txt
	python ./savage_dickey.py 0.01 "M_Cocle->DiegoAguacate" 1 10 ./run$i/mcmc7.txt
	python ./savage_dickey.py 0.01 "M_DiegoAguacate->Cocle" 1 10 ./run$i/mcmc7.txt
	python ./savage_dickey.py 0.01 "M_Piedras->DiegoAguacate" 1 10 ./run$i/mcmc8.txt
	python ./savage_dickey.py 0.01 "M_DiegoAguacate->Piedras" 1 10 ./run$i/mcmc8.txt
	python ./savage_dickey.py 0.01 "M_Piedras->DiegoAguacateCocle" 1 10 ./run$i/mcmc9.txt
	python ./savage_dickey.py 0.01 "M_DiegoAguacateCocle->Piedras" 1 10 ./run$i/mcmc9.txt
done