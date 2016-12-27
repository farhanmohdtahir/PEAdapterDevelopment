#include "CS.h"
#include <string.h>
#include <cstring>
#include <string>
#include <iostream>
        
using namespace std;			

const char CS::nuclist[4] = {'A','C','G','T'};

CS::CS()
{
	memset(nucleotidecount,0,sizeof(nucleotidecount));
	memset(phred,0,sizeof(phred));
	memset(consensus,0,sizeof(consensus));
	adapterLength = 0;
	adapterPos = 0;
	addPrior = 0;
}

void CS::cs(string seq_1, int max, int & c1)      
{ 						                     	               	
        c1 = 0;
        
        if(addPrior == 0)
	{
		for(int i = 0; i < 20; i++)
		{
			for(int j = 0; j < 4; j++)
			{
				nucleotidecount[j][i] += 0.25;
			}
		}
		addPrior++;
	}
	
        for(int c = max; c < seq_1.length(); c++)
        {
	   if(c1 == 20) break;
           switch(seq_1[c])
           {
               case 'A': nucleotidecount [0][c1] += 1;
                         break;
               case 'C': nucleotidecount [1][c1] += 1;
                         break;
               case 'G': nucleotidecount [2][c1] += 1;
                         break;
               case 'T': nucleotidecount [3][c1] += 1;
                         break;
	       case 'N': nucleotidecount [0][c1] += 1;
			 nucleotidecount [1][c1] += 1;
			 nucleotidecount [2][c1] += 1;
			 nucleotidecount [3][c1] += 1;
           }
	   c1++;
        }
        
        for(int d = 0; d < 20; d++)
        {
           int max = -1, count = 0;
           for(int e = 0; e < 4; e++)
           {   if ( e == 0) max = 0;
               if(nucleotidecount[e][d] > nucleotidecount[max][d])
               { 
                   max = e;

               } 
           }

           if(max > -1)
           {
               count++;
           }
	
	 if(count == 1)                                             /*Determining nucleotide IUPAC */
           {                                                        /*symbol for consensus sequence.*/ 
		if(max == 0)consensus[d] = 'A';
		if(max == 1)consensus[d] = 'C';
		if(max == 2)consensus[d] = 'G';		
		if(max == 3)consensus[d] = 'T';
           }
	 else  
           {
               consensus[d] = 'N';
           }
        }
        
	calc_phred();
}

void CS::calc_phred()
{	
	for(int i = 0; i < 20; i++)
	{
		float total1 = (nucleotidecount[0][i]+nucleotidecount[1][i]+nucleotidecount[2][i]+nucleotidecount[3][i]);
		if(nucleotidecount[0][i] != 0)  
			phred[0][i] = (-10*log10(1-(nucleotidecount[0][i]/total1)));
		if(nucleotidecount[1][i] != 0)
			phred[1][i] = (-10*log10(1-(nucleotidecount[1][i]/total1)));
		if(nucleotidecount[2][i] != 0) 
			phred[2][i] = (-10*log10(1-(nucleotidecount[2][i]/total1)));
		if(nucleotidecount[3][i] != 0)
			phred[3][i] = (-10*log10(1-(nucleotidecount[3][i]/total1)));		
	}
}

void CS::print_nucCount_phred()
{
        cout << "\nNucleotide Count :\n"; 

        cout << endl;
        for(int p = 0; p < 4; p++)
        {
                cout << nuclist[p] << " ";
                for(int q = 0; q < 20; q++)
                { 
                        if(q > 8) cout << " ";
                                cout << setw(5) << nucleotidecount [p][q] << " ";
                }
                cout << endl;
        }
        cout << endl;
        
        cout << "\nPhred score :\n"; 

        cout << endl;
        for(int p = 0; p < 4; p++)
        {
                cout << nuclist[p] << " ";
                for(int q = 0; q < 20; q++)
                { 
                        if(q > 8) cout << " ";
                        	cout.precision(4);
                                cout << setw(9) << phred [p][q] << " ";
                }
                cout << endl;
        }
        cout << endl;
}

void CS::print_cs(int c, int opt = 0)
{
        if(opt == 1)
        {
        	for(int complement = 0; complement < c; complement++)
        	{

              		 if(consensus[complement] == 'A')
             		         consensus[complement] = 'T';
            		 else if(consensus[complement] == 'T')
             		         consensus[complement] = 'A';
             		 else if(consensus[complement] == 'C')
             		         consensus[complement] = 'G';
              		 else if(consensus[complement] == 'G')
                 		 consensus[complement] = 'C';
      		}
        }

        for(int r = 0; r < c; r++)
        {
        	bool noZero = false;	
        	for(int s = 0; s < 4; s++)
        	{
        		if(nucleotidecount[s][r] != 0)
        			noZero = true;
        	}
        	if(noZero == true)
                	cout << consensus[r]; 
        }
}

void CS::checkConfidence(double conf, int& confTrue, int adapLenCount)
{
	int checkLength = 0;
int checkPos = 0;

	for(int a = 0; a < 20; a++)
	{
		Confidence[a]=0;
		bool highConf = false;
		double biggestValue = 0;
		for(int b = 0; b < 4; b++)
		{
			if(phred[b][a] != 0)
			{
			if(phred[b][a] > biggestValue)
					biggestValue = phred[b][a];
			}
			else checkLength += 1; 
		}
		if(biggestValue >= conf)
			Confidence[a] = 1; 
		
		if(checkLength == 4)
		{
			checkPos = a;
		}
		checkLength = 0;
	}
	
	if(adapLenCount > 0)
	{
		if(adapterPos == checkPos)
		adapterLength++;
		else 
		{
			adapterPos = checkPos;	
			adapterLength = 0;
		}
	}
	else adapterPos = checkPos;
	
	if(adapterLength >= 5)
	{	
		if(adapterPos == 0)
			adapterPos = 20;
		int count = 0;  
		for(int c = 0; c < adapterPos; c++)
		{
			if(Confidence[c] == 1)
				count++;		
		}
	
		if(count == adapterPos)
		{
			confTrue++;
		}
		adapterPos = 0;
	}
	
} 