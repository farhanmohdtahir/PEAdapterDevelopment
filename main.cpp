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
    
    
	string file1 = "", file2 = "", seq_1, seq_2, seq_2_rev, seq_1_al, seq_2_al, dash1, dash2;
  	int opt = 0, seqLength = 0, debugLevel = 0;
  	double percentage = .0, confLevel = .0;
  	bool option = false;
  	int rowmax = 0, colmax = 0, confTrue, adapLenCount = 0, iteration = 0;
  	bool skip = true, fourline1, fourline2;
        int bil=1, count=0, idline=0, dnaline=1, c1, c2, adapLen1=0, adapLen2=0;
        
 	string line, line2;
 	bool onlynuc = true;
 	bool onlynuc2 = true;
 	
 	  static struct option long_options[] = {
        {"help",                     no_argument,       0,  'h' },
        {"f1",                   required_argument,     0,  'a' },
        {"f2",                   required_argument,     0,  'b' },
        {"lengthPercentage",     required_argument,     0,  'l' },
        {"matchPercentage",      required_argument,     0,  'm' }, 
        {"confidenceLevel",      required_argument,     0,  'c' },  
	{"debugLevel",           required_argument,     0,  'd' },                                
        {0,                               0,            0,   0  }
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
             case 'l' : seqLength = atoi(optarg);option = true;
                 break;
             case 'm' : percentage = atoi(optarg);option = true;
                 break;  
             case 'c' : confLevel = atoi(optarg);option = true;
                 break; 
             case 'd' : debugLevel = atoi(optarg);option = true;
                 break;                                                      
             default: print_usage(); 
                 exit(EXIT_FAILURE);
        }
    }

    if(seqLength == 0) seqLength = 70;
    if(percentage == 0) percentage = 85;
    if(confLevel == 0) confLevel = 10;
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

                    ab.complementInput(seq_2);
                    seq_2_rev=seq_2;                                        //assign reverse-complement read 2 into variable seq_2_rev
                    reverse( seq_2_rev.begin(), seq_2_rev.end() );          //assign reverse-complement read 2 into variable seq_2_rev
                    

                    double L1 = seq_1.length();
                    double L2 = seq_2_rev.length();	   

                    b.nw(seq_1, seq_2_rev, seq_1_al, seq_2_al, debugLevel);

                    if(b.percentage > percentage && (b.rowmax/L1*100) > seqLength && ((seq_2_rev.length()-b.colmax)/L2*100) > seqLength && b.colmax != 0)
                    {
                        
                        rowmax = b.rowmax;
                        colmax = seq_2_rev.length()-b.colmax;                   
                        c.cs(seq_1, rowmax, c1);
                        d.cs(seq_2, colmax, c2); 
                        
                        if(debugLevel == 1 || debugLevel == 2){
                            dash1="";
                            dash2="";
                            
                            for(int i=0; i<c1; i++){
                                dash1+="-";
                            }
                            for(int i=0; i<c2; i++){
                                dash2+="-";
                            }                   
                            cout<<"\nSeq1:"<<dash2<<seq_1<<endl;
                            cout<<"Seq2:"<<seq_2_rev<<dash1<<endl;
                        }
                        
                            
                        if(debugLevel == 1 || debugLevel == 2) c.print_nucCount_phred();
                        if(debugLevel == 1 || debugLevel == 2) d.print_nucCount_phred();
                        
                        confTrue = 0;
                        c.checkConfidence(confLevel, confTrue, adapLenCount, adapLen1);
                        d.checkConfidence(confLevel, confTrue, adapLenCount, adapLen2);
                        
                        adapLenCount++;                        
                        if(confTrue == 2)
                        {
                                cout <<endl<<"Adapter in FASTQ file 1: ";
                                c.print_cs(adapLen1, 0);
                                cout <<endl<<"Adapter in FASTQ file 2: ";
                                d.print_cs(adapLen2, 1);
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
    cout<<"Usage: ./PEAdapterFinder -f1 filename1 -f2 filename2 -l lengthPercentage -m matchPercentage -c confidenceLevel -d debugLevel\n";
}

void help()
{
    cout<<"\nPaired-End Adapter Finder\n\n";
    cout<<"Please enter ./PEAdapterFinder -f1 filename1 -f2 filename2 -l lengthPercentage -m matchPercentage -c confidenceLevel -d debugLevel\n";
    cout<<"\n-f1 [required argument] - name of first fastq file \n";
    cout<<"-f2 [required argument] - name of second fastq file\n";
    cout<<"-l  [optional argument] - minimum length percentage to get adapter sequence (default = 70)\n";
    cout<<"-m  [optional argument] - minimum match percentage to get adapter sequence (default = 85)\n";   
    cout<<"-c  [optional argument] - minimum confidence level of nucleotides (default = 10)\n";
    cout<<"-d  [optional argument] - debug level of programme (default = 0 : 0 - only adapter sequences, 1 - nucleotide count and phred score, 3 - dynamic programming matrix and traceback matrix)\n";         
    cout<<"\nPlease refer to the documentation for more help...\n";
    cout<<"\nThank you. \n\n";
}