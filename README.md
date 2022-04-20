# Explicit Caching HYB based SpMV
This is the codes for paper: [Explicit caching HYB: a new high-performance SpMV framework on GPGPU](https://arxiv.org/abs/2204.06666)
The code is tested in the matrices listed in the Appendix of that paper

Please compile the code according to your environment configurations. Change the $CUDA\_HOME variable in the Makefile.
In my experiment, I built a directory named "read" at the project directory, and store the .mtx files in that directory, if you want to store the .mtx files in other path, please modify the "fopen" function arguments at the solver\_test.c file. 

After the code is compiled, users can test the code as follow:
```
./spmv.out -i 2000 -m audikw_1
```

Where -i indicates the number of iterations for test, and the -m indicates the file name of the matrix for test (the matrix should be stored with a .mtx file)
