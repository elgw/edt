CC=cc -std=c99
CFLAGS=-Wall -Wextra -O3 -pedantic -fopenmp
LDFLAGS=-flto -lm -fopenmp

FILES=eudist.o src/eudist_cli.c
eudist: $(FILES)
	$(CC) $(CFLAGS) $(FILES) $(LDFLAGS) -o eudist

eudist.o: src/eudist.c
	$(CC) -c $(CFLAGS) src/eudist.c

clean:
	rm -f eudist
	rm -f eudist.o

python:
	python3 eudist_setup.py build_ext --inplace
	python3 eudist_setup.py install
