#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#include <cstdio>
#include <cmath>
#include "NistTests.h"	//for macros

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       M A T R I X   R A N K  A L G O R I T H M  F U N C T I O N  
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define	MATRIX_FORWARD_ELIMINATION	0
#define	MATRIX_BACKWARD_ELIMINATION	1

static void perform_elementary_row_operations(int flag, int i, int M, int Q, char **A)
{
	int		 j;
	int		 k;
	
	if ( flag == MATRIX_FORWARD_ELIMINATION ) {
		for ( j=i+1; j<M;  j++ )
			if ( A[j][i] == 1 ) 
				for ( k=i; k<Q; k++ ) 
					A[j][k] = (A[j][k] + A[i][k]) % 2;
	}
	else {
		for ( j=i-1; j>=0;  j-- )
			if ( A[j][i] == 1 )
				for ( k=0; k<Q; k++ )
					A[j][k] = (A[j][k] + A[i][k]) % 2;
	}
}

static int swap_rows(int i, int index, int Q, char **A)
{
	int		p;
    char	temp;
	
	for ( p=0; p<Q; p++ ) {
		temp = A[i][p];
		A[i][p] = A[index][p];
		A[index][p] = temp;
	}
	
	return 1;
}

static int find_unit_element_and_swap(int flag, int i, int M, int Q, char **A)
{ 
	int		 index;
	int		 row_op=0;
	
	if ( flag == MATRIX_FORWARD_ELIMINATION ) {
		index = i+1;
		while ( (index < M) && (A[index][i] == 0) ) 
			index++;
			if ( index < M )
				row_op = swap_rows(i, index, Q, A);
	}
	else {
		index = i-1;
		while ( (index >= 0) && (A[index][i] == 0) ) 
			index--;
			if ( index >= 0 )
				row_op = swap_rows(i, index, Q, A);
	}
	
	return row_op;
}

static int determine_rank(int m, int M, int Q, char **A)
{
	int		 i;
	int		 j;
	int		 rank;
	int		 allZeroes;
	
	/* DETERMINE RANK, THAT IS, COUNT THE NUMBER OF NONZERO ROWS */
	
	rank = m;
	for ( i=0; i<M; i++ ) {
		allZeroes = 1; 
		for ( j=0; j<Q; j++)  {
			if ( A[i][j] == 1 ) {
				allZeroes = 0;
				break;
			}
		}
		if ( allZeroes == 1 )
			rank--;
	} 
	
	return rank;
}

int computeRank(int M, int Q, char **matrix)
{
	int		 i;
	int		 rank;
	int		 m=MIN(M,Q);
	
	/* FORWARD APPLICATION OF ELEMENTARY ROW OPERATIONS */ 
	for ( i=0; i<m-1; i++ ) {
		if ( matrix[i][i] == 1 ) 
			perform_elementary_row_operations(MATRIX_FORWARD_ELIMINATION, i, M, Q, matrix);
		else { 	/* matrix[i][i] = 0 */
			if ( find_unit_element_and_swap(MATRIX_FORWARD_ELIMINATION, i, M, Q, matrix) == 1 ) 
				perform_elementary_row_operations(MATRIX_FORWARD_ELIMINATION, i, M, Q, matrix);
		}
	}

	/* BACKWARD APPLICATION OF ELEMENTARY ROW OPERATIONS */ 
	for ( i=m-1; i>0; i-- ) {
		if ( matrix[i][i] == 1 )
			perform_elementary_row_operations(MATRIX_BACKWARD_ELIMINATION, i, M, Q, matrix);
		else { 	/* matrix[i][i] = 0 */
			if ( find_unit_element_and_swap(MATRIX_BACKWARD_ELIMINATION, i, M, Q, matrix) == 1 )
				perform_elementary_row_operations(MATRIX_BACKWARD_ELIMINATION, i, M, Q, matrix);
		}
	} 

	rank = determine_rank(m, M, Q, matrix);

	return rank;
}

