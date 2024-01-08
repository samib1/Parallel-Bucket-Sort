# Parallelized Bucket Sort

**Purpose:** To implement a shared memory and distributed computing version of parallel bucket sort

**Author:** Sami

**Program Files included:**  README.md, shared_mem.cpp, distributed.cpp

---


# Overview
Bucket sort is used to sort non-uniform distributed numbers by placing them into buckets of equal number of points per bucket regardless of the distribution of values. 

Here accelerated the serial bucket sort using a distributed and shared memory approach utilizing OpenMPI and OpenMP respectively.

---

# How To Run The Program
## 1. Running OpenMPI (for distributed)
1. compile with ``` mpicxx  -std=c++17  -o output distributed.cpp ```

2. Run command => ```mpirun -np num_processors -oversubscribe output input_numbers num_samples```.
Sample run command => ```mpirun  -np 4 -oversubscribe output 100000 10 ``` in cmdline.

3. Clean up with ```rm output```


## 2. Running OpenMP (for shared memory)
1. compile with ``` g++  -std=c++17  -fopenmp  -o output shared_mem.cpp```

2. Run command => ```./output input_numbers num_samples num_threads>>```.
Sample run command => ```./output 100000 10 6 >>./part3 num samples threads``` in cmdline.

3. Clean up with ```rm output```

## 3. Ooutput
1. The time to sort the samples will be printed on the commandline

---
## NOTES
1. To run this program you must have c++, OpenMPI, openMP installed.

___

**Have a wonderful day!**
