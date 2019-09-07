test.exe: testperm.f90 ../src/signtestperm.f90
	gfortran -fbounds-check testperm.f90 ../src/signtestperm.f90 -o test.exe
clean:
	/bin/rm -f *.o *.exe ../src/*.o ../src/*.so
veryclean:
	find . -name "*.o" -delete -o -name "*.so" -delete -o -name ".R*" -delete
probest.so: ../src/probest.f90
	R CMD SHLIB ../src/probest.f90
