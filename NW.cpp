#include "NW.h"
#include <string.h>
#include "cmath"

using namespace std;



NW::NW()
{
	F = NULL;
	Fx = 0;
	Fy = 0;
	rowmax = 0;	
	colmax = 0;	
	count = 0;
}

NW::~NW()
{
	//clear();
}

int NW::nw(string seq_1, string seq_2, string& seq_1_al, string& seq_2_al, int debug)
{
        int  d = 2 ;                 /* Gap penalty */

	seq_1_al = "";
	seq_2_al = "";

        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

if(count == 0)
{ 
   	F = new int * [ L2+1 ];
       	for( int i = 0; i <= L2; i++ )  F[ i ] = new int [ L1 ];

      	// Traceback matrix
       	traceback = new char * [ L2+1 ];
       	for( int i = 0; i <= L2; i++ )  traceback[ i ] = new char [ L1 ];

       	        
        count++;
}

if(L1 > Fx || L2 > Fy)
{       
	//clear();
	// Dynamic programming matrix
   	F = new int * [ L2+1 ];
       	for( int i = 0; i <= L2; i++ )  F[ i ] = new int [ L1 ];

      	// Traceback matrix
       	traceback = new char * [ L2+1 ];
       	for( int i = 0; i <= L2; i++ )  traceback[ i ] = new char [ L1 ];
}

	Fx = L1;
	Fy = L2;
	
        // Initialize traceback and F matrix (fill in first row and column)
        dpm_init( F, traceback, L1, L2, d );

        // Create alignment
        nw_align( F, traceback, seq_1, seq_2, seq_1_al, seq_2_al, d);
	
	if(debug == 2)
	{
        	cout << "\nDynamic Programming Matrix: " << "\n\n";
        	print_matrix( F, seq_1, seq_2 );

        	cout << "\nTraceback Matrix: " << "\n\n";
        	print_traceback( traceback, seq_1, seq_2, d, L1, L2 );

		cout << "\nAligned Sequences: " << "\n\n";
        	print_al( seq_1_al, seq_2_al );  
        }
        
        return  0 ;
}


void  NW::dpm_init( int ** F, char ** traceback, int L1, int L2, int d )
{
        F[ 0 ][ 0 ] =  0 ;
        traceback[ 0 ][ 0 ] = 'n' ;

        int i=0, j=0;

        for( j = 1; j <= L1; j++ )
        {
                F[ 0 ][ j ] =   -j * d ;
                traceback[ 0 ][ j ] =  '-' ;
        }
        for( i = 1; i <= L2; i++ )
        {
                F[ i ][ 0 ] =  0;
                traceback[ i ][ 0 ] =  '|' ;
        }
}

int NW::nw_align(                   
               int **     F,
              char **    traceback,
              string     seq_1,
              string     seq_2,
              string&    seq_1_al,
              string&    seq_2_al,
              int        d
            )
{
        int         L1 = seq_1.length();
        int         L2 = seq_2.length(); 
        int checkMatch=0; 
        int        fU;
        int        fD;
        int        fL;
        char       ptr;
        char       nuc[L2], nuc2[L1] ;
        int        i = 0, j = 0, k = -1*10^9999;

        const int  a =  2;   /* Match */
        const int  b = -1;   /* Mismatch */


        
        strncpy(nuc, seq_1.c_str(), sizeof(nuc));
        strncpy(nuc2, seq_2.c_str(), sizeof(nuc2));
        
        for( i = 0; i < L2; i++ )
        {
                for( j = 0; j < L1; j++ )
                {  
                    if (nuc[j]!='N'||nuc2[i]!='N'){    //Determine match or mismatch between 2 bases
                        if (nuc[j]==nuc2[i]){
                            checkMatch=a;
                        }
                        else {
                            checkMatch=b;
                        }
                    }
                        else {
                        checkMatch=0; 
                        }
                        
                        fU = F[ i ][ j+1 ] - d ;
                        fD = F[ i ][ j ] + checkMatch ;
                        fL = F[ i+1 ][ j ] - d ;
                        
                        F[ i+1 ][ j+1 ] = max(max(fU, fD),fL);  //finding the highest score between fU, fD and fL

                }
        }
        
        for(int m = 0; m <= L1; m++)
        {
                if( F[ L2 ][ m ] > k )
                {
                        k = F[ L2 ][ m ];
                        rowmax = m;
                }
        }
        j = rowmax;

//store traceback matrix
        while( i > 0  && j > 0 )
        {
            if (nuc[j-1]!='N'||nuc2[i-1]!='N'){    
                        if (nuc[j-1]==nuc2[i-1]){
                            checkMatch=a;
                        }
                        else {
                            checkMatch=b;
                        }
                    }
                        else {
                        checkMatch=0; 
                        }
            
            if (F[i][j]==F[i-1][j]-d){
                seq_1_al += '-' ;            
                seq_2_al += seq_2[ i-1 ] ; 
                i-- ;
            }
            
            else if (F[i][j]==F[i-1][j-1]+checkMatch){
                seq_1_al += seq_1[ j-1 ] ; 
                seq_2_al += seq_2[ i-1 ] ; 
                i-- ;  j-- ; 
            }
            
            else {
                seq_1_al += seq_1[ j-1 ] ; 
                seq_2_al += '-' ; 
                j--;
            }
	
		if(j==0 )
		{
			colmax = i;
			break;						
		}
            
            if(i==0){
                break;
            }
        }


	verifyPercentage(seq_1_al, seq_2_al );
	
        reverse( seq_1_al.begin(), seq_1_al.end() );
        reverse( seq_2_al.begin(), seq_2_al.end() );

        return  0 ;
}

void  NW::print_matrix( int ** F, string seq_1, string seq_2 )
{
        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

        cout << "        ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << "   ";
        }
        cout << "\n  ";

        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout.width( 3 );
                        cout << F[ i ][ j ] << " ";
                }
                cout << endl;
        }
}

void  NW::print_traceback( char ** traceback, string seq_1, string seq_2, int d, int L1, int L2 )
{
    int checkMatch;
    char nuc[L1], nuc2[L2] ;
    const int  a = 2, b= -1;   
    
    strncpy(nuc, seq_1.c_str(), sizeof(nuc));
    strncpy(nuc2, seq_2.c_str(), sizeof(nuc2));
        
        cout << "    ";
        
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << " ";
        }
        cout << endl<< "  ";
        
        for( int i = 0; i <=L2 ; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                
                for(int j=0; j<=L1; j++){

                    if (i==0 || j==0){
                        cout<<traceback[i][j]<<" ";
                    }
                    
                    else{
                    
                    if (nuc[j-1]!='N'||nuc2[i-1]!='N'){    
                            if (nuc[j-1]==nuc2[i-1]){
                                checkMatch=a;
                            }
                            else {
                                checkMatch=b;
                            }
                        }
                            else {
                            checkMatch=0; 
                            }

                    if (F[i][j]==F[i-1][j-1]+checkMatch){
                        cout<<'\\'<<" ";
                    }
                    
                    else if (F[i][j]==F[i-1][j]-d){
                        cout<<'|'<< " ";
                    }


                    else {
                        cout<<'-'<<" ";
                    }
                }
                    
              }
              cout<<endl;
            }
}

void  NW::print_al( string& seq_1_al, string& seq_2_al )
{
        cout << seq_1_al << endl;
        cout << seq_2_al << endl;
}

void  NW::verifyPercentage( string seq_1_al, string seq_2_al )
{
	int AL1 = seq_1_al.length();
	double score = .0;
	percentage = .0;
	for(int i = 0; i < AL1; i++)
	{
		if(seq_1_al[i]==seq_2_al[i])
			score++;
	}
	percentage = (score/AL1)*100;	
}

void NW::clear()
{
	for( int i = 0; i <= Fy; i++ )  delete[] F[ i ];
        delete[] F;
        for( int i = 0; i <= Fy; i++ )  delete[] traceback[ i ];
        delete[] traceback;
}
