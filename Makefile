eudist:
	gcc -Wall -Wextra -pedantic-errors -O3 -std=c11 -march=native eudist.c -lm -flto -lpthread -o eudist

test:
	gcc -Wall -std=c99 -g eudist.c -lm -lpthread -o eudist_test
	valgrind -v ./eudist_test
#-fopenmp

python:
	python3 eudist_setup.py build_ext --inplace
	python3 eudist_setup.py install
