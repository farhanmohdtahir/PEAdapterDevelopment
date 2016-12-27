
#include "iostream"
#include "string"
#include "algorithm"
#include "stdlib.h"
#include "stdio.h"
//#include "emmintrin.h"
//#include "xmmintrin.h"



#ifndef NW_H
#define NW_H

/** @brief The NW Class corresponds to the Needleman-Wunsch algorithm which is used to align the two reads.  

	The Needleman-Wunsch algorithm is implemented thorought dynamic programming. Read 1 is aligned against the reverse complement of Read 2 to obtain the best alignment. The nucleotides at the right end of Read 1 and upper end of Read 2 form the adapter sequences for Adapter 1 and Adapter 2 respectively.
	
	@author Rayan Gan
	@date April 2015    
 */
using namespace std;
class NW
{
private:
	/** Initialize the dynamic programming matrix and traceback matrix
	 * @param F The dynamic programming matrix
	 * @param traceback The traceback matrix
	 * @param L1 The length of Read 1
	 * @param L2 The length of Read 2
	 * @param d Gap penalty
	 */
	void  dpm_init( int ** F, char ** traceback, int L1, int L2, int d );
	/** Runs the Needleman-Wunsch Alignment to get the best alignment
	 * @param F The dynamic programming matrix
	 * @param traceback The traceback matrix
	 * @param seq_1 Read 1 in the first fastq file
	 * @param seq_2 Read 2 in the second fastq file
	 * @param seq_1_al Read 1 aligned to Read 2
	 * @param seq_2_al Read 2 aligned to Read 1
	 * @param d Gap penalty
	 */
	int nw_align(int ** F, char ** traceback, string seq_1, string seq_2, string& seq_1_al, string& seq_2_al, int d);
	/** Prints out the dynamic programming matrix
	 * @param F The dynamic programming matrix
	 * @param seq_1 Read 1 in the first fastq file
	 * @param seq_2 Read 2 in the second fastq file
	 */
	void  print_matrix( int ** F, string seq_1, string seq_2 );
        /** Prints out the traceback matrix
	 * @param traceback The traceback matrix
	 * @param seq_1 Read 1 in the first fastq file
	 * @param seq_2 Read 2 in the second fastq file	
	 */
	void  print_traceback( char ** traceback, string seq_1, string seq_2, int d, int L1, int L2);
	/** Prints out the aligned sequences
	 * @param seq_1_al Read 1 aligned to Read 2
	 * @param seq_2_al Read 2 aligned to Read 1	
	 */
	void  print_al( string& seq_1_al, string&  seq_2_al );
	/** Determines the traceback percentage
	 * @param seq_1_al Read 1 aligned to Read 2
	 * @param seq_2_al Read 2 aligned to Read 1	
	 */
	void verifyPercentage( string seq_1_al, string seq_2_al );
	
	/** Clears the dynamic programming and traceback matrices
	 */
	void clear();
        
public:
	/**Main function in the NW class which calls other functions to perform dynamic programming on Reads 1 and Reads 2
	 * @param seq_1 Read 1 in the first fastq file
	 * @param seq_2 Read 2 in the second fastq file
	 * @param seq_1_al Read 1 aligned to Read 2
	 * @param seq_2_al Read 2 aligned to Read 1
	 */
	int nw(string seq_1, string seq_2, string& seq_1_al, string& seq_2_al, int debug);	
	
	int rowmax;		/**< Largest value at last row of dynamic programming. */
	int colmax;		/**< Traceback value at first column of dynamic programming. */
	int count;		/**< Counts the number of runs of NW. */
	int Fx;			/**< Keeps track of the length of sequence 1. */
	int Fy;			/**< Keeps track of the length of sequence 2. */
	int ** F;		/**< Stores the dynamic programming values. */
	char ** traceback;	/**< Stores the traceback characters. */
	char * tracebackscore;	/**< Stores the best traceback characters. */
	double percentage;	/**< Percentage of matches in best traceback. */

	NW();
	~NW();
};



#endif