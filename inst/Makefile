allexec: testconcor.exe testperm.exe testaov.exe testtskmsurv.exe
	touch testaov.exe
testconcor.exe: testconcord.f90 ../src/ufcn.f90
	gfortran -Wall -fbounds-check -fcheck=all testconcord.f90 ../src/ufcn.f90 -o testconcord.exe
testperm.exe: testperm.f90 ../src/signtestperm.f90
	gfortran -Wall -fbounds-check -fcheck=all testperm.f90 ../src/signtestperm.f90 -o test.exe
testaov.exe: testaov.o aovp.o nextp.o
	gfortran -Wall -fbounds-check -fcheck=all -o testaov.exe testaov.o aovp.o nextp.o
testtskmsurv.exe: tskmsurv.o testtskmsurv.o rankem.o tsksurvpv.o nextp.o ../src/initperm.f90
	gfortran -Wall -fbounds-check -fcheck=all -o testtskmsurv.exe tskmsurv.o testtskmsurv.o rankem.o tsksurvpv.o nextp.o ../src/initperm.f90
nextp.o: ../src/nextp.f90
	gfortran -Wall -fbounds-check -fcheck=all -c ../src/nextp.f90
aovp.o: ../src/aovp.f90
	gfortran -Wall -fbounds-check -fcheck=all -c ../src/aovp.f90
clean:
	/bin/rm -f *.o *.exe ../src/*.o ../src/*.so
veryclean:
	find . -name "*.o" -delete -o -name "*.so" -delete -o -name ".R*" -delete
probest.so: ../src/probest.f90
	R CMD SHLIB ../src/probest.f90

tsksurvpv.o: ../src/tsksurvpv.f90
	gfortran -Wall -fbounds-check -fcheck=all -c ../src/tsksurvpv.f90
testtskmsurv.o: testtskmsurv.f90
	gfortran -Wall -fbounds-check -fcheck=all -c -o testtskmsurv.o testtskmsurv.f90
rankem.o: ../src/rankem.f90
	gfortran -Wall -fbounds-check -fcheck=all -c -o rankem.o ../src/rankem.f90
#tskmsurv.o: ../src/tskmsurv.f90
#	gfortran -Wall -fbounds-check -fcheck=all -c -o tskmsurv.o ../src/tskmsurv.f90
tskmsurv.o: tskmsurv.F90
	gfortran -Wall -fbounds-check -fcheck=all -c -o tskmsurv.o tskmsurv.F90
testaov.o: testaov.f90
	gfortran -Wall -fbounds-check -fcheck=all -c -o testaov.o testaov.f90
