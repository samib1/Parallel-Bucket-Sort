/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Author:     Sami Byaruhanga 
Purpose:    To implement OpenMP shared mem solution to sample sort
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


//INCLUDES + GLOBAL VARS ********************************************
#include <iostream>
#include <omp.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>
#include <unistd.h>

using namespace std;
//************************************************* INCLUDES END HERE


//FUNCTIONS  ********************************************************
//-------------------------------------------------------
// printVector
// PURPOSE: Prints the vector information
// INPUT PARAMETERS: vector and its name
//------------------------------------------
void printVector(std::vector<int> inputVec, string vecName){
    //-------------------------------------
    std::cout << vecName <<" => ";
    for (int num : inputVec) {
        std::cout << num << " ";
    }
    std::cout << std::endl;
    //-------------------------------------
}
//-------------------------------------------------------
//************************************************ FUNCTIONS END HERE


//MAIN **************************************************************
int main(int argc, char**argv){
    //Initializations ----------------------------------
    if (argc < 4) {
        std::cerr << "Provide an array length argument and number of samples" << std::endl;
        std::exit(-1);
    }
    int N = atoi(argv[1]);
    int nSamples = atoi(argv[2]); // int nSamples = 10; //can change this 
    int nthreads = atoi(argv[3]);
    int nBuckets = nthreads;
    int nSplitters = nBuckets -1; 
    if (nSamples < nSplitters){
        std::cerr << "Number of samples, s, should not be less than number of splitters, p." << std::endl;
        exit(1);
    }
    //--------------------------------------------------

    //Generate my list ---------------------------------
    std::random_device rd{};
    std::mt19937 generator{rd()};
    std::poisson_distribution distribution(4); //If an event occurs 4 times a minute on average, how often is it that it occurs n times in one minute?
    // std::normal_distribution distribution{500.0, 100.0};
    // std::binomial_distribution distribution(4,0.5); //perform 4 trials, each succeeds 1 in 2 times
    // std::geometric_distribution distribution; // Kinda like coin tosses
    //function definition for generating random numbers
    auto random_int = [&distribution, &generator]{ return std::round(distribution(generator)); };    
    std::vector<int> values(N);
    for (int n = 0; n < N; n++){
        values[n] = random_int();
    }
    //--------------------------------------------------

    
    //Timing -------------------------------------------
    omp_set_num_threads(nthreads);
    double tstartParallel = omp_get_wtime();
    //--------------------------------------------------
    
    
//??LIKELY PUT THIS IN PARALLEL FOR AND HAVE SINGLES AND ALL FOR??
    //1. multithread: sample ---------------------------
    //generate samples by sampling the list uniformly
    std::default_random_engine generator2;
    std::uniform_int_distribution<int> sample_distribution(0, N-1);
    std::vector<int> samples(nSamples);
    #pragma omp parallel for 
    for (int isample = 0; isample < nSamples; isample++)
    {
        samples[isample] = *(values.begin() + sample_distribution(generator2));
    }
    // printVector(samples, "Samples");
    //--------------------------------------------------
    

    //2. Single: sort and pick p-1 splitters -----------
    //put the splitters in order so they can be used to define buckets 
    sort(samples.begin(), samples.end());
// ??DO WE NEED A FOR HERE?? WELL NEED ONLY SINGLE TO DO THIS    
    std::vector<int> splitters(nSplitters);
    for (int isplitter = 0; isplitter < nBuckets; isplitter++)
    {
        splitters[isplitter] = samples[(isplitter+1)*nSamples / nBuckets];
    }
    // printVector(splitters, "Splitters");
    //--------------------------------------------------


    //3. Multithread: partion and sort -----------------
    std::vector<std::vector<int>::iterator> bucket_bounds(nBuckets + 1);  
    bucket_bounds[0] = values.begin();  
    for (int isplitter = 0; isplitter < nSplitters; isplitter++)
    {
        //partition function is useful here
        bucket_bounds[isplitter+1] = std::partition(bucket_bounds[isplitter], values.end(), [&](int value) { return value < splitters[isplitter]; } );
    }
    bucket_bounds[nBuckets] = values.end();
    
    //sort buckets
    #pragma omp parallel for
    for (int ibucket = 0; ibucket < nBuckets; ibucket++)
    {
        sort(bucket_bounds[ibucket], bucket_bounds[ibucket+1]);

        int id = omp_get_thread_num();        
        // sleep(id + 1); //just checking if it prints 
        // std::cout << "id " << id << " of " << nthreads << std::endl;
    }
    //--------------------------------------------------
    

    //--------------------------------------------------
    double tstopParallel = omp_get_wtime();
    std::cout << "Parallel Time = " << tstopParallel - tstartParallel << " seconds." << std::endl;
    //--------------------------------------------------

    //printing list ------------------------------------
    // printVector(values, "Full unsorted list");
    // std::cout << "\nSorted list => ";
    // std::copy(values.begin(), values.end(), std::ostream_iterator<int>(std::cout, " "));
    // std::cout << " " <<std::endl;
    return 0;    
}
//**************************************************** MAIN ENDS HERE


/* STEPS TO SOLVE PROBLEM *******************************************
 1. multithread: sample
 2. Single: sort and pick p-1 splitters 
 3. Multithread: partion and sort
*********************************************************************/
