function eudist_build()
mex ../src/eudist_mex.c ../src/eudist.c ...
    -output eudist CFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
end