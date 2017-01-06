#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <getopt.h>
#include "Input.h"
#include "NW.h"
#include "CS.h"
#include "time.h"

using namespace std;

/**	@mainpage Paired-End Adapter Finder v15.4.14

	The Paired-End Adapter Finder constructs two consensus adapter sequences from Reads 1 and Reads 2 of paired-end sequencing by using the Needleman-Wunsh algorithm to align the two reads and determining which region of the sequence is actually an adapter. The user has to input two fastq format files, corresponding to Read 1 and Read 2, and the result will be the consensus adapter sequences for both the Adapters in Read 1 and Read 2 respectively of the paired-end sequencing.
	
	@author Rayan Gan and Farhan Tahir
	@date January 2017
 */
 
void print_usage();
void help();

int main(int argc, char *argv[]) 
{
    
	string file1 = "", file2 = "", seq_1, seq_2, seq_1_al, seq_2_al;
  	int opt = 0, seqLength = 0, debugLevel = 0;
  	double percentage = .0, confLevel = .0;
  	bool option = false;
  	int rowmax = 0, colmax = 0, confTrue, adapLenCount = 0, iteration = 0;
  	bool skip = true, fourline1, fourline2;
        int bil=1, count=0, idline=0, dnaline=1, adapMax1=0, adapMax2=0, c1, c2;
        
 	string line, line2;
 	bool onlynuc = true;
 	bool onlynuc2 = true;
 	
 	  static struct option long_options[] = {
        {"help",      no_argument,       0,  'h' },
        {"f1",    required_argument,     0,  'a' },
        {"f2",    required_argument,     0,  'b' },
        {"seql",  optional_argument,     0,  'c' },
        {"perc",  optional_argument,     0,  'd' }, 
        {"conf",  optional_argument,     0,  'e' },  
	{"debug", optional_argument,     0,  'f' },                                
        {0,                0,            0,   0  }
    };

    int long_index =0;
    while ((opt = getopt_long_only(argc, argv,"", 
                   long_options, &long_index )) != -1) {
        switch (opt) {
             case 'h' : help();option = true;
                 break;
             case 'a' : file1 = optarg;option = true;
                 break;
             case 'b' : file2 = optarg;option = true;
                 break;
             case 'c' : seqLength = atoi(optarg);option = true;
                 break;
             case 'd' : percentage = atoi(optarg);option = true;
                 break;  
             case 'e' : confLevel = atoi(optarg);option = true;
                 break; 
             case 'f' : debugLevel = atoi(optarg);option = true;
                 break;                                                      
             default: print_usage(); 
                 exit(EXIT_FAILURE);
        }
    }
    
    if(seqLength == 0) seqLength = 70;
    if(percentage == 0) percentage = 85;
    if(confLevel == 0) confLevel = 1;
    if(debugLevel == 0) debugLevel = 0;

    if(option == false) print_usage();
	if(debugLevel == 0 || debugLevel == 1 || debugLevel == 2){
    if(file1 != "" && file2 != "")
    {
    	Input ab;
	NW b;
 	CS c;
	CS d;
 	ifstream myfile (file1.c_str());
 	ifstream myfile2 (file2.c_str());
        ofstream outfile1;
        ofstream outfile2; 

//  Determine whether both Input file is multi-line FASTQ file or 4-line FASTQ file
 	if (myfile.is_open() && myfile2.is_open())
  	{
            while (getline (myfile,line)) 
            {
                   if (count==idline){  
                       if (line[0]=='@'){
                           fourline1=true;
                           idline+=4;
                       }
                       else{
                           fourline1=false;
                           break;
                       }
                   }
               ++count;
               if (count>=200) break;
            }        

         count=0; idline=0;
         
         while(getline(myfile2,line2)){
                if (count==idline){  
                    if (line2[0]=='@'){
                        fourline2=true;
                        idline+=4;
                    }
                    else{
                        fourline2=false;
                        break;
                    }
                }
            ++count;
            if (count>=200) break;
         }
         myfile.close();
         myfile2.close();
        }
         
        else {
            cout << "Unable to open file"<<endl; 
            exit(0);
        }         
         
        if (fourline1==true&&fourline2==true) cout<<"This is normal 4-line FASTQ file. "<<endl<<endl;
          
        else cout<<"This is multi-line FASTQ file."<<endl<<"Reformation into 4-Line FASTQ file..."<<endl<<endl;
         
//Reformation of multi-line FASTQ file into 4-line FASTQ file     
        if (fourline1==false) file1=ab.reform(file1, fourline1);

        if (fourline2==false) file2=ab.reform(file2, fourline2);
 
//Finding Adapter sequence
        if (fourline1==true && fourline2==true){ 
        ifstream infile1(file1.c_str());
        ifstream infile2(file2.c_str());
        count=1;
        
        cout<<"Finding adapter..."<<endl;
        while (getline (infile1,line) && getline (infile2,line2))
        {
                     if(count==dnaline){
                                
                        infile1>>seq_1;
                        infile2>>seq_2;
                        dnaline+=4;
                        
// Creating NW and CS objects                    
                    reverse( seq_2.begin(), seq_2.end() );
                    ab.complementInput(seq_2);

                    double L1 = seq_1.length();
                    double L2 = seq_2.length();	   

                    b.nw(seq_1, seq_2, seq_1_al, seq_2_al, debugLevel);

                    if(b.percentage > percentage && (b.rowmax/L1*100) > seqLength && ((seq_2.length()-b.colmax)/L2*100) > seqLength && b.colmax != 0)
                    {
                        if(debugLevel == 1 || debugLevel == 2)
                        {
                            cout<<"\nSeq1:"<<seq_1<<endl;
                            cout<<"Seq2:"<<seq_2<<endl;
                        }

                        rowmax = b.rowmax;
                        colmax = seq_2.length()-b.colmax;                    
                        c.cs(seq_1, rowmax, c1);
                        reverse( seq_2.begin(), seq_2.end() );
                        d.cs(seq_2, colmax, c2); 

                        if(debugLevel == 1 || debugLevel == 2)c.print_nucCount_phred();;
                        if(debugLevel == 1 || debugLevel == 2)d.print_nucCount_phred();

                        confTrue = 0;
                        c.checkConfidence(confLevel, confTrue, adapLenCount);
                        d.checkConfidence(confLevel, confTrue, adapLenCount);
                        adapLenCount++;                        
                        if (c1>adapMax1) adapMax1=c1;
                        if (c2>adapMax2) adapMax2=c2;

                        if(confTrue == 2)
                        {
                                cout <<endl<<"Adapter in FASTQ file 1: ";
                                c.print_cs(adapMax1, 0);
                                cout <<endl<<"Adapter in FASTQ file 2: ";
                                d.print_cs(adapMax2, 1);
                                cout<<endl;
                                exit(0);
                        }	
                    }

                    ++bil;
        // After NW and CS
                    }

                    ++count;
                 }
  	  infile1.close();
  	  infile2.close();
  	  cout << "Confidence level could not be achieved...\n"; 
        }
    }
}
    
  return 0;
}

void print_usage() 
{
    cout<<"Usage: ./PEAdapterFinder -f1 filename -f2 filename -seql percentage length -perc percentage match -conf confidence level -debug debug level\n";
}

void help()
{
    cout<<"\nPaired-End Adapter Finder\n\n";
    cout<<"Please enter ./PEAdapterFinder -f1 file1 -f2 file2 -seql seqLength -perc percentage -conf confLevel -debug debugLevel\n";
    cout<<"\nfile1 - name of first fastq file \n";
    cout<<"file2 - name of second fastq file\n";
    cout<<"seqLength - minimum length percentage to get adapter sequence (default = 70, to change use '-seql=')\n";
    cout<<"percentage - minimum match percentage to get adapter sequence (default = 85, to change use '-perc=')\n";   
    cout<<"confLevel - minimum confidence level of nucleotides (default = 1, to change use '-conf=')\n";
    cout<<"debugLevel - debug level of programme (default = 0, to change use '-debug=' : 0 - only adapter sequences, 1 - nucleotide count and phred score, 3 - dynamic programming matrix and traceback matrix)\n";         
    cout<<"\nPlease refer to the documentation for more help...\n";
    cout<<"\nThank you. \n\n";
}