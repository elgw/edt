eudist:
	gcc -Wall -O3 -std=c11 -march=native eudist.c -lm -flto  -fopenmp -DOMP

eudist_no_openmp:
	gcc -Wall -O3 -std=c11 -march=native eudist.c -lm -flto 


test:
	gcc -Wall -std=c99 -g eudist.c -lm -o eudist_test
	valgrind -v ./eudist_test
#-fopenmp

