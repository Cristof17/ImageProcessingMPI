out_file=filtru
topo=topology.txt
imgs=imagini.txt
out=statistica.txt

build:
	mpicc echo.c -o $(out_file)

run:
	./$(out_file) $(topo) $(imgs) $(out) 

clean:
	rm -rf 343C1_RotschingCristofor_Tema3.zip
	rm -rf filtru

zip:
	zip 343C1_RotschingCristofor_Tema3.zip *.c *.h makefile Readme

