
#!/bin/sh

mpiexecjl --project=../ -n 4 julia -J ../Channel.so -O3 --check-bounds=no -e 'using Channel; Channel.main(16,1,4)'
