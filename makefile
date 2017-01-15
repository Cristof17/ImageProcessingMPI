out_file=filtru
topo=topology.txt
imgs=imagini_test.txt
out=statistica.txt
processes=12

build:
	mpicc echo.c -o $(out_file)

run:
	mpirun -np $(processes) $(out_file) $(topo) $(imgs) $(out) 

clean:
	rm -rf 343C1_RotschingCristofor_Tema3.zip
	rm -rf filtru

zip:
	zip 343C1_RotschingCristofor_Tema3.zip *.c *.h makefile Readme

