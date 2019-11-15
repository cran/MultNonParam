testconcord.exe: testconcord.f90 ../src/ufcn.f90
	gfortran -fbounds-check testconcord.f90 ../src/ufcn.f90 -o testconcord.exe
testperm.exe: testperm.f90 ../src/signtestperm.f90
	gfortran -fbounds-check testperm.f90 ../src/signtestperm.f90 -o test.exe
testaov.exe: testaov.f90 aovp.o nextp.o
	gfortran -fbounds-check -o testaov.exe testaov.f90 aovp.o nextp.o
nextp.o: ../src/nextp.f90
	gfortran -fbounds-check -c ../src/nextp.f90
aovp.o: ../src/aovp.f90
	gfortran -fbounds-check -c ../src/aovp.f90
clean:
	/bin/rm -f *.o *.exe ../src/*.o ../src/*.so
veryclean:
	find . -name "*.o" -delete -o -name "*.so" -delete -o -name ".R*" -delete
probest.so: ../src/probest.f90
	R CMD SHLIB ../src/probest.f90
