function eudist_build()
mex('../src/eudist_mex.c', '../src/eudist.c', ...
    '-output' , 'eudist')
end