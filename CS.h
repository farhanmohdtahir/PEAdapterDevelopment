#include "math.h"
#include "iomanip"
#include "string"

#ifndef CS_H
#define CS_H

/** @brief The CS Class corresponds to the Consensus Sequence which is used obtain the adapter sequences from the two reads.      

	The Consensus Sequence for both Read 1 and Read 2 adapters are determined based on the traceback of the dynamic programming algorithm. The Consensus Sequence is searched based on how well Read 1 and Read 2 align to each other. The ends of Read 1 and Read 2 which do not align form adapters and are added to the nucleotide count. The final Consensus Sequence is based on the nucleotide count.
	
	@author Rayan Gan
	@date April 2015
 */
using namespace std;
class CS
{
private:
	static const char nuclist[4];		/**< List of nucleotides. */	

	double nucleotidecount [4][20];		/**< 2D-array of number of nucleotides at each posiiton in the first adapter. */
	double phred [4][20];			/**< 2D-array of Phred score of each nucleotide in the first adapter. */
	char consensus[20]; 			/**< Consensus sequence of first adapter. */
	int Confidence[20];
	int adapterLength;
	int adapterPos;	
	int addPrior;

public:
	/**Main function in CS class which performs the algorithm to find the consensus sequences 
	 * @param seq_1 Read 1 in the first fastq file
	 * @param max Length at which alignment is the best	 
	 */
	void cs(string seq_1, int max, int & c1 );
	/**Calculates the Phred scores of each nucleotide in the nucleotide count
	 */
	void calc_phred();
	/**Prints out the Consensus Sequence
	 * @param opt Option to choose if Read 1 or Read 2 (0 - Read 1, 1 - Read 2)
	 */
	void print_cs(int c, int opt);
	/**Prints out the nucleotide count and Phred scores
	 */	
	void print_nucCount_phred();
	/**Checks the confidence level of the consensus sequences 
	 * @param conf Confidence level
	 * @param confTrue Does the sequence meet the confidence level?	 
	 */
	void checkConfidence(double conf, int& confTrue, int adapLenCount);
	CS(); 
};

#endif