test.exe: testperm.f90 ../src/signtestperm.f90
	gfortran -fbounds-check testperm.f90 ../src/signtestperm.f90 -o test.exe
clean:
	/bin/rm -f *.o *.exe ../src/*.o ../src/*.so
