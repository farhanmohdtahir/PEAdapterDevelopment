#include "iostream"
#include "fstream"
#include "string"
#include "algorithm"
#include "stdlib.h"
#include "stdio.h"
#include "Input.h"

using namespace std;


/** @brief The Input Class transforms the sequences to be used in the Adapter Sequencer.      

	The Input Class is responsible to reform multi-line FASTQ file into 4-line FASTQ file (in case) and reverse complement Read 2 before dynamic programming is done. 
	@author Rayan Gan and Farhan Tahir
	@date January 2017
 */


/**Reverse complement the input for sequences from the second input file.     
 * @param seq The sequence to be complemented
 */
int Input::complementInput(string& seq)
{
        for(int reverse = 0; reverse < seq.length(); reverse++)
        {
            switch(seq[reverse]){
                case 'A': seq[reverse]='T';break;
                case 'T': seq[reverse]='A';break;
                case 'G': seq[reverse]='C';break;
                case 'C': seq[reverse]='G';
            }
        }      
	return 0;
}

//Function to reformat the multi-line FASTQ file into 4-line FASTQ file
//@param file Input File, fourline Determine whether the multi-line FASTQ file completely transform into fourLine FASTQ file

string Input::reform(string file, bool &fourline){
    bool startline=false;
    ifstream myfile;
    ofstream outfile;
    string line, str;
    int cmpvar=0;
    
    myfile.open(file.c_str(),ios::in);
    file="4line_"+file;
    outfile.open(file.c_str(),ios::out); 
    
    while (getline(myfile, line)){
        
                if (startline==true) {
                        if(line[0]=='@'){
                            if (cmpvar==0) str+=line;
                            if (cmpvar==2) {
                            str+="\n"+line;
                            cmpvar=3;
                            }           
                            else{
                            str+="\n"+line;
                            cmpvar=0;
                            }
                        }     
                        
                        else if(line[0]=='A'||line[0]=='T'||line[0]=='G'||line[0]=='C'||line[0]=='N'){                           
                            if (cmpvar==1||cmpvar==3) str+=line;
                            else if (cmpvar==0) {
                                str+="\n"+line;
                                cmpvar=1;
                            }
                            else if (cmpvar==2) {
                                str+="\n"+line;
                                cmpvar=3;
                            }
                        }
                        
                        else if(line[0]=='+'){
                            if (cmpvar==2) {
                                str+="\n"+line;
                                cmpvar=3;
                            }
                            else{
                            str+="\n"+line;
                            cmpvar=2;
                            }
                        }
                        
                        else {
                            if (cmpvar==0) str+=line;
                            
                            else if (cmpvar==2) {
                                str+="\n"+line;
                                cmpvar=3;
                            }
                            else str+=line;
                            
                        }
                    }
                    
                else {
                        if(line[0]=='@'){
                            startline=true;
                            str+=line;
                            cmpvar=0;
                        }
                    }
    }
            outfile<<str;
            outfile.close();
            fourline=true;
            return file;
}

