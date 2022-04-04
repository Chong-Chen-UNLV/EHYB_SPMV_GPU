# Explicit Caching HYB based SpMV
This is the codes for paper: Explicit caching HYB: a new high-performance SpMV framework on GPGPU 
The code is tested in the matrices listed in the Appendix of that paper

to test the file please fist derived the cache memory size using equation (1) and (2) of that paper, 
for example for matrix audikw_1, the dimension of matrix is 77651847, we can derive the size 
of the vector cache in the shared memory: 93KB

The we change the line 23 of the kernel.h (memPerThread)to 93, that is, make the line 93 of kernel.h
```
const int memPerThread = 93;
```
generate the executable file:
```
make
```

Then we can test the file as follow:
```
./spmv.out -i 2000 -m audikw_1
```

Where -i indicates the number of iterations for test, and the -m indicates the matrix for test
where matrix should stored in the directory: ../read as a .mtx file (there should be a directory
named read at the fathor directory of the executable file, and the matrix file should be a 
.mtx file)
