/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Author:     Sami Byaruhanga 
Purpose:    Implement distributed version of sample sort
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


//INCLUDES + GLOBAL VARS ********************************************
#include <iostream>
#include <mpi.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>
using namespace std;
//************************************************* INCLUDES END HERE


//FUNCTIONS  ********************************************************
// parallelRange ----------------------------------------
// PURPOSE: Tells how data partioned among p/sors
// INPUT PARAMETERS: Ordered list of seq entries w/ 1st, last element
// RETURNS: indexes (Local start, stop) and
//          t.t num of local entries to be processed
//------------------------------------------
void parallel_range(int rank, int nproc, int globalStart, int globalStop, int &localStart, int &localStop, int &localCount){
	int globalSize = globalStop - globalStart;
	int localSize = globalSize / nproc;
	int localRemainder = globalSize % nproc;

    int offset = rank < localRemainder ? rank : localRemainder;
	
	localStart = rank*localSize + offset;
	localStop = localStart + localSize;
	if(localRemainder > rank) localStop++;

	localStart += globalStart;
	localStop += globalStart;
    localCount = localStop - localStart;
}
//-------------------------------------------------------


// printVector ------------------------------------------
// PURPOSE: Prints the vector information
// INPUT PARAMETERS: vector and its name
//------------------------------------------
void printVector(std::vector<int> inputVec, string vecName){
    //-------------------------------------
    std::cout << vecName <<" => ";
    for (int num : inputVec) {
        std::cout << num << " ";
    }
    std::cout << "\n" << std::endl;
    //-------------------------------------
}
//-------------------------------------------------------
//************************************************ FUNCTIONS END HERE


//MAIN **************************************************************
int main(int argc, char**argv){
    //Initializations ---------------------------------------------
    int rank, nproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
    if (argc < 3) {
        std::cerr << "Provide an array length argument and number of samples" << std::endl;
        std::exit(-1);
    }
    int N = atoi(argv[1]);
    int nSamples = atoi(argv[2]); // int nSamples = 10; //can change this 
    // int N = std::stoi(argv[1]); //array size to sort 
    // int nSamples = std::atoi(argv[2]); // int nSamples = 10; //can change this 
    int nBuckets = nproc;
    int nSplitters = nBuckets -1; 
    if (nSamples < nSplitters){
        std::cerr << "Number of samples, s, should not be less than number of splitters, p." << std::endl;
        exit(1);
    }
    //-------------------------------------------------------------

    if(rank ==0){
        //0. Random number generation from a distribution ---------
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
        // printVector(values, "full unSorted");
        
        //0.1 Distrubition so can scatter to all ranks 
        std::vector<int> localStartBuffer(nproc), localCountBuffer(nproc), localStopBuffer(nproc);
        for(int i=0; i<nproc;i++){
            parallel_range(i, nproc, 0, N, localStartBuffer[i], localStopBuffer[i], 
            localCountBuffer[i]);
        }
        std::vector<int> unsorted(localCountBuffer[0]);
        MPI_Scatterv(&values[0], &localCountBuffer[0], &localStartBuffer[0], MPI_INT, &unsorted[0], 
            localCountBuffer[0], MPI_INT, 0, MPI_COMM_WORLD); 
        // std::string nameArr = "Rank " + std::to_string(rank) + " array";
        // printVector(unsorted, nameArr);
        //---------------------------------------------------------


        //START TIMMING -------------------------------------------
        MPI_Barrier(MPI_COMM_WORLD);
        double tStart = MPI_Wtime();
        //----------------------------------------------------------

        //1. Sample the values on each procesor [use PROF CODE] ---
        //generate samples by sampling the list uniformly
        std::default_random_engine generator1;
        std::uniform_int_distribution<int> sample_distribution(0, localCountBuffer[0]-1);
        std::vector<int> samples(nSamples);
        for (int isample = 0; isample < nSamples; isample++){
            samples[isample] = *(unsorted.begin() + sample_distribution(generator1));
        }
        // std::string nameSamplesArr = "Samples on rank " + std::to_string(rank);
        // printVector(samples, nameSamplesArr);
        //---------------------------------------------------------

        //@ rank 0 gather sample ----------------------------------
        // MPI_Barrier(MPI_COMM_WORLD);
        std::vector<int> allSamples(nSamples*nproc); //all ranks have same num of samples so get allsamples size
        MPI_Gather(&samples[0], nSamples, MPI_INT, &allSamples[0], nSamples, MPI_INT, 0, MPI_COMM_WORLD);
        // printVector(allSamples, "All samples on rank 0");

        //2. Sort samples 
        sort(allSamples.begin(), allSamples.end());
        // printVector(allSamples, "All samples sorted on rank 0");

        //3. Pick p-1 spitters & broadcast splitters
        std::vector<int> splitters(nSplitters); // use nbuckets - 1 for splitters
        for (int isplitter = 0; isplitter < nSplitters; isplitter++){
            splitters[isplitter] = allSamples[(isplitter+1)*nSamples*nproc / nBuckets];
        }
        printVector(splitters, "splitters rank 0");        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&splitters[0], nSplitters, MPI_INT, 0, MPI_COMM_WORLD);
        //-------------------------------------------------------------


        //4. All processors -------------------------------------------
        //Put w/ partitions in buckets each p/sor & distribute [ALTOALLV]
        // printVector(unsorted, "RANK 0 BEFORE PARTION");
        std::vector<std::vector<int>::iterator> bucket_bounds(nBuckets + 1);  
        bucket_bounds[0] = unsorted.begin();  
        for (int isplitter = 0; isplitter < nSplitters; isplitter++){
            bucket_bounds[isplitter+1] = std::partition(bucket_bounds[isplitter], unsorted.end(), [&](int value) { return value < splitters[isplitter]; } );
        }
        bucket_bounds[nBuckets] = unsorted.end();
        // printVector(unsorted, "RANK 0 AFTER PARTION"); //IF IT WORKS WE SHUD HAVE SPLITTERS AS BOUNDS e.g., if its 2,5,8 splitters, our list should have vales less than 2 on one side and should not apper after 5 and so on
        
        //Now distribute to all processors LUCAS CODE AS GUIDE
        // setup alltoallv for unsorted to sorted
        // MPI_Barrier(MPI_COMM_WORLD);
        std::vector<int> sendcounts(nproc);
        std::vector<int> senddispls(nproc);
        std::vector<int> recvcounts(nproc);
        std::vector<int> recvdispls(nproc);

        //setup sendcounts
        for(int ibucket=0; ibucket < std::size(bucket_bounds); ibucket++) {
            sendcounts[ibucket] = std::distance(bucket_bounds[ibucket], bucket_bounds[ibucket+1]);
        }

        // alltoall turns sendcounts into recvcounts
        MPI_Alltoall(&sendcounts[0], 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);

        // setup senddispls and recvdispls
        // displs[0] are zero because we want values to start at
        // sorted[0] on ever process
        senddispls[0] = 0;
        recvdispls[0] = 0;

        for(int ibucket=1; ibucket < std::size(bucket_bounds); ibucket++) {
            senddispls[ibucket] = senddispls[ibucket - 1] + sendcounts[ibucket - 1];
            recvdispls[ibucket] = recvdispls[ibucket - 1] + recvcounts[ibucket - 1];
        }

        // 3) last step is to resize sorted so that it can contain all received values
//??MIGHT NEED TO JUST REUSE UNSORTED BUT FOR NOW WE CREATE A SORTED VECTOR
        std::vector<int> sorted;
        sorted.resize(recvcounts[nproc-1] + recvdispls[nproc-1]);
        // std::cout << "Rank " << rank << " sending " << unsorted.size() << " values. Receiving " << sorted.size() << " values." << std::endl;
        std::string be4AllV = "Be4 AllV Rank " + std::to_string(rank);
        // printVector(unsorted, be4AllV);

        MPI_Alltoallv(&unsorted[0], &sendcounts[0], &senddispls[0], MPI_INT,
                        &sorted[0], &recvcounts[0], &recvdispls[0], MPI_INT, MPI_COMM_WORLD);
        std::string afterAllV2 = "After AllV Rank " + std::to_string(rank);
        // printVector(sorted, afterAllV2);
        //-------------------------------------------------------------


        //5. Sort buckets responsible for  ----------------------------
        sort(sorted.begin(), sorted.end());
        std::string afterSorting = "Sorted on Rank " + std::to_string(rank) + " ";
        // printVector(sorted, afterSorting);
        //-------------------------------------------------------------

        //Timing information ------------------------------------------
        MPI_Barrier(MPI_COMM_WORLD);
        double tStop = MPI_Wtime();
        std::cout << "Time = " << tStop - tStart << " seconds" << std::endl;
        //-------------------------------------------------------------
    }
    else{
        //0. Receiving random list ------------------------------------
        int localStart, localStop, localCount;
        parallel_range(rank, nproc, 0, N, localStart, localStop, localCount);
        std::vector<int> unsorted(localCount);
        MPI_Scatterv(NULL, NULL, NULL, NULL, &unsorted[0], localCount, MPI_INT, 0, MPI_COMM_WORLD);
        std::string nameArr = "Rank " + std::to_string(rank) + " array";
        // printVector(unsorted, nameArr);
        //-------------------------------------------------------------

        //SYNC for timing ---------------------------------------------
        MPI_Barrier(MPI_COMM_WORLD);
        //-------------------------------------------------------------

        //1. Sample the values on each procesor [use PROF CODE] -------
        //generate samples by sampling the list uniformly
        std::default_random_engine generator1;
        std::uniform_int_distribution<int> sample_distribution(0, localCount-1);
        std::vector<int> samples(nSamples);
        for (int isample = 0; isample < nSamples; isample++){
            samples[isample] = *(unsorted.begin() + sample_distribution(generator1));
        }
        
        //send all samples to rank 0
        // MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(&samples[0], nSamples, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
        //-------------------------------------------------------------


        //All processors ----------------------------------------------
        MPI_Barrier(MPI_COMM_WORLD);
        std::vector<int> splitters(nSplitters); // use nbuckets - 1 for splitters
        MPI_Bcast(&splitters[0], nSplitters, MPI_INT, 0, MPI_COMM_WORLD);
        // printVector(splitters, "SPITTER RECVED");

        //4. Put with partitions in buckets each p/sor and distribute [ALTOALLV]
        // std::string be4Partion = "Rank " + std::to_string(rank) + " be4 partition array";
        // printVector(unsorted, be4Partion);
        std::vector<std::vector<int>::iterator> bucket_bounds(nBuckets + 1);  
        bucket_bounds[0] = unsorted.begin();  
        for (int isplitter = 0; isplitter < nSplitters; isplitter++){
            //partition function is useful here
            bucket_bounds[isplitter+1] = std::partition(bucket_bounds[isplitter], unsorted.end(), [&](int value) { return value < splitters[isplitter]; } );
        }
        bucket_bounds[nBuckets] = unsorted.end();
        // std::string afterPartion = "Rank " + std::to_string(rank) + " after partition array";
        // printVector(unsorted, afterPartion);

        //Now distribute to all processors 
        // MPI_Barrier(MPI_COMM_WORLD);
        std::vector<int> sendcounts(nproc);
        std::vector<int> senddispls(nproc);
        std::vector<int> recvcounts(nproc);
        std::vector<int> recvdispls(nproc);

        //setup sendcounts
        for(int ibucket=0; ibucket < std::size(bucket_bounds); ibucket++) {
            sendcounts[ibucket] = std::distance(bucket_bounds[ibucket], bucket_bounds[ibucket+1]);
        }
        // printVector(sendcounts, "Send counts");

        // alltoall turns sendcounts into recvcounts
        MPI_Alltoall(&sendcounts[0], 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);

        // setup senddispls and recvdispls
        senddispls[0] = 0;
        recvdispls[0] = 0;

        for(int ibucket=1; ibucket < std::size(bucket_bounds); ibucket++) {
            senddispls[ibucket] = senddispls[ibucket - 1] + sendcounts[ibucket - 1];
            recvdispls[ibucket] = recvdispls[ibucket - 1] + recvcounts[ibucket - 1];
        }

        // 3) last step is to resize sorted so that it can contain all received values
//??MIGHT NEED TO JUST REUSE UNSORTED BUT FOR NOW WE CREATE A SORTED VECTOR
        std::vector<int> sorted;
        sorted.resize(recvcounts[nproc-1] + recvdispls[nproc-1]);
        // std::cout << "Rank " << rank << " sending " << unsorted.size() << " values. Receiving " << sorted.size() << " values." << std::endl;

        std::string be4AllV = "Be4 AllV Rank " + std::to_string(rank);
        // printVector(unsorted, be4AllV);
        MPI_Alltoallv(&unsorted[0], &sendcounts[0], &senddispls[0], MPI_INT,
                        &sorted[0], &recvcounts[0], &recvdispls[0], MPI_INT, MPI_COMM_WORLD);
        std::string afterAllV2 = "After AllV Rank " + std::to_string(rank);
        // printVector(sorted, afterAllV2);
        //-------------------------------------------------------------

        //5. Sort buckets responsible for  ----------------------------
        sort(sorted.begin(), sorted.end());
        std::string afterSorting = "Sorted on Rank " + std::to_string(rank) + " ";
        // printVector(sorted, afterSorting);
        //-------------------------------------------------------------

        //Timing ------------------------------------------------------
        MPI_Barrier(MPI_COMM_WORLD);
        //-------------------------------------------------------------
    }    
    MPI_Finalize();
}
//**************************************************** MAIN ENDS HERE


/** NOTES ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 * QUESTIONS: 
 * 1.ON SPLITTERS: wont the buckets just uneven distributions???
     ON STEP3: Pick p-1 spitters //???IF I BROADCAST AND USE THAT TO PARTITION WILL END UP HAVING SOME WITH NO VALUSE
 *                  
 * 
 * IMRPOVEMENTS:
 * 1. AM DOING TWICE THE WORK SO DOING ONLY THE RANK0 WORK IN IF AND ALL ELSE OUTSIDE will make code clean
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/


/** STEPS TO SOLVE PROGRAM *******************************************
 * 0. Generate random #s [Non-uniform] on rank 0 and distribute
 * 1. Sample the values on each procesor [PROF CODE]
 * @ rank 0 gather
 *      2. Sort samples 
 *      3. Pick p-1 spitters
 *      broadcast
 * All processors
 * 4. Put with partitions in buckets each p/sor and distribute [ALTOALLV]
 * 5. Sort buckets responsible for  
*********************************************************************/
