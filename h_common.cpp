#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>

#include "constants.h"
#include "externs.h"
#include "h_externs.h"
#include "common.h"
#include "h_struct.h"
#include "h_common.h"
#include "params.h"
#include "simfold.h"

// Hosna feb 12, 2008
#include "W_final.h"
#include "hfold.h"
// Hosna, May 3rd, 2012
//#include "hfold_pkonly.h"
// Hosna, November 16, 2015
#include "hfold_interacting.h"

//kevin
#include <vector>
#include <utility>
#include "init.h"

/*
 * This function is just the same as detect_original_pairs
 * but it also calculates the arcs for each base and saves them in arc_table
 *
 * The algorithm for finding the arcs is as follows:
 * When checking each base in the input structure,
 *
 * case 1) if (structure[i] == '.' or ' ' or '_')
 * 		1-a) if stack is not empty, then put arc_table[i] = top element on the stack
 * 		1-b) otherwise arc_table[i] = -1
 *
 * case 2) if (structure[i] == '(')
 * 		2-a) if stack is not empty, then put arc_table[i] = top element on the stack and push i in the stack
 * 		2-b) otherwise arc_table[i] = -1 and push i in the stack
 *
 * case 3) if (strcuture[i] == ')'), pop the stack and match it with j
 * 		3-a) if stack is not empty, then put arc_table[i] = top element on the stack
 * 		3-b) otherwise put arc_table[i] = -1
 *
 */

void detect_original_pairs_arcs(char *structure, int *p_table, int *arc_table)
// PRE:  structure contains the desired structure
// POST: pairs will contain the index of each base pair
//               or -1 if it does not pair
//
{
		//printf("structure: %s\n", structure);
        int i, j, struct_len;
        stack_ds st;
        h_init (&st);
        remove_space (structure);
		//printf("The given structure is: \n %s\n", structure);
        struct_len = strlen (structure);
	// Hosna March 8, 2012
	// since index i starts at 0 and stack top is also set to 0 to show stack is empty, if i=0 is paired then we have incorrect arc values!
	// So I am introducing STACK_EMPTY = -1 to h_common.h and change h_init and h_pop accordingly

        for (i=0; i < struct_len; i++)
          {
			  // Hosna March 8, 2012
			  // changing nested ifs to switch for optimality
			  switch (structure[i])
				{
					case '.':
					{
					  p_table[i] = RESTRICTED_UNPAIR;
					  if (st.top > STACK_EMPTY){//0){
		//              	if (debug)
		//	              	printf("base %d: stack has %d elements in it and its top element is %d \n",i, st.top, st.elem[st.top]);
						arc_table[i] = st.elem[st.top];
					  }else{
						arc_table[i] = NOT_COVERED;
					  }
					}
						break;
					case ' ':
					case '_':
					{
					  p_table[i] = FREE_TO_PAIR;
					  if (st.top > STACK_EMPTY){//0){
		//              	if (debug)
		//	              	printf("base %d: stack has %d elements in it and its top element is %d \n",i, st.top, st.elem[st.top]);
						arc_table[i] = st.elem[st.top];
					  }else{
						arc_table[i] = NOT_COVERED;
					  }
					}
						break;
					case '(':
						{
							if (st.top > STACK_EMPTY){//0){
			//              	if (debug)
			//	              	printf("base %d: stack has %d elements in it and its top element is %d \n",i, st.top, st.elem[st.top]);
								arc_table[i] = st.elem[st.top];
							  }else{
								arc_table[i] = NOT_COVERED;
							  }
							  h_push (&st, i);
						}
						break;
					case ')':
				//else if (structure[i] == ')')
					  {
						j = h_pop (&st);
						p_table[i] = j;
						p_table[j] = i;
						  if (st.top > STACK_EMPTY){//0){
		//                	if (debug)
		//	                	printf("base %d: stack has %d elements in it and its top element is %d \n",i, st.top, st.elem[st.top]);
							arc_table[i] = st.elem[st.top];
							}else{
								arc_table[i] = NOT_COVERED;
							}
					  }
						break;
			  }
          }

	/* for (i=137; i<struct_len; i++){
			printf("p_table[%d] = %d AND arc_table[%d]=%d \n\n",i,p_table[i],i,arc_table[i]);
		}
		printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"); */
        if (st.top != STACK_EMPTY)//0)
        {
            fprintf (stderr, "The given structure is not valid: %d more left parentheses than right parentheses\n", st.top);
            exit (1);
        }

}

/**
 * This function is just like the above function except the case that
 * it can handle the pseudoknotted structures of at most density 2
 *
 * the density 2 structures can be presented with ( and [ in dot parenthesis format
 * so we need at most 2 stacks to keep track of the base pairings
 * we call the stacks st and st_brack for ( and [ respectively.
 *
 *
 */

void detect_original_PKed_pairs(char *structure, int *p_table)
// PRE:  structure contains the desired structure
// POST: pairs will contain the index of each base pair
//               or -1 if it does not pair
{
        int i, j, struct_len;
        stack_ds st; //stach used for (
        stack_ds st_brack; // stack used for [
        h_init (&st);
        h_init (&st_brack);
        remove_space (structure);
        struct_len = strlen (structure);
        for (i=0; i < struct_len; i++)
          {
			  // Hosna March 8, 2012
			  // changing nested ifs to switch for optimality
			  switch (structure[i])
			  {
				  case '.':
            //if (structure[i] == '.')
					{
					  p_table[i] = -1;
					}
					  break;
				  case ' ':
				  case '_':
            //else if (structure[i] == ' ' || structure[i] == '_')
					{
					  p_table[i] = -2;
					}
					  break;
				  case '(':
            //else if (structure[i] == '(')
					{
					  h_push (&st, i);
					}
					  break;
				  case '[':
            //else if (structure[i] == '[')
					{
					  h_push (&st_brack, i);
					}
					  break;
				  case ')':
            //else if (structure[i] == ')')
					{
						j = h_pop (&st);
						p_table[i] = j;
						p_table[j] = i;
					}
					  break;
				  case ']':
           // else if (structure[i] == ']')
					{
						j = h_pop (&st_brack);
						p_table[i] = j;
						p_table[j] = i;
					}
					  break;
			  }
          }
        if (st.top != STACK_EMPTY || st_brack.top != STACK_EMPTY) //0 || st_brack.top != 0)
        {
            fprintf (stderr, "The given structure is not valid: %d more left parenthesis than right parentheses\n", st.top);
            exit (1);
        }
}

/*
 * The algorithm for finding the weakly closed regions for a given pseudoknot free structure is as follows
 *
 * initialization:
 * 				open = -1;
 * 				weakly_closed array initialized to 0
 *
 *
 * case 1) if the region is of the form [i,i] and i is not paired,
 * then it is considered as a weakly closed region and we put weakly_closed[ii] = 1
 * otherwise we put open = i
 *
 * case 2) if we are considering region [i,j] where i!=j and we know region [i, j-1] is weakly closed
 *		2-a) if j is not paired, then [i,j] is also weakly closed and we put weakly_closed[ij] = 1
 * 		2-b) if j is paired and we have i <= bp(j) < j then [i,j] is weakly_closed and we put weakly_closed[ij] = 1
 * 		THIS CASE CAN NEVER HAPPERN ==> REMOVE FROM CODE
 * 		2-c) otherwise [i,j] is NOT weakly closed
 * 			if (open == -1) then open = j
 *
 * case 3) if we are considering region [i,j] where i!=j and we know region [i,j-1] is NOT weakly closed
 * 		3-a) if j is paired and bp(j) == open then [i,j] is weakly closed and we put weakly_closed[ij] = 1 and open = -1
 * 		3-b) otherwise the region is still NOT weakly closed
 *
 */

void detect_weakly_closed(h_str_features *fres, int *weakly_closed, int nb_nucleotides, int *index){
	int i,j;
	int open = -1;

	for (i = 0; i < nb_nucleotides; i++){
		open = -1;
		for (j = i; j < nb_nucleotides; j++){
			int ij = index[i]+j-i; // index[i]+j-i gives the index ij
//			if ( i == 7 && j == 27){
//				printf("weakly_closed[%d,%d] = %d and open = %d \n", i,j-1,weakly_closed[ij-1], open);
//			}
			if (i == j ){
				if (fres[j].pair < 0){ // j is not paired
					weakly_closed[ij]= 1;
//					if (debug ){//_weakly_closed){
//						printf("%d==%d and %d is not paired => weakly closed\n",i,j,j);
//					}
				}else{
					open = j;
//					if (debug){//_weakly_closed){
//						printf("%d==%d and %d is paired => not weakly closed, and open = %d\n",i,j,j,open);
//					}
				}
			}else if (weakly_closed[ij-1]){
				if (fres[j].pair < 0){ // j is not paired
					weakly_closed[ij] = 1;
//					if (debug ){//_weakly_closed){
//						printf("[%d,%d] is weakly closed and %d is not paired => weakly closed\n",i,j-1,j);
//					}
				}else if (open == -1){
					open = j;
//					if (debug ){//_weakly_closed){
//						printf("[%d,%d] is weakly closed but %d is paired => not weakly closed, and open = %d\n",i,j-1,j,open);
//					}
				}
			}else if(fres[j].pair >= 0 && open == fres[j].pair){
				weakly_closed[ij] = 1;
				open = -1;
//				if (debug && i == 27 && j == 68){//_weakly_closed){
//					printf("[%d,%d] is not weakly closed but %d is paired with %d => weakly closed\n",i,j-1,j,fres[j].pair);
//				}
			}
			//if (debug){//_weakly_closed){
//				if (i == 27 && j == 68){//_weakly_closed){
//					printf("weakly_closed[%d,%d] = %d\n",i,j,weakly_closed[ij]);
//				}
//				if (weakly_closed[ij] == 1)
//					printf("Region [%d,%d] is Weakly closed\n",i,j);
			//}


		}
	}
}
/*
 * Hosna: Feb 12, 2007
 * algorithm for finding if region [i,j] is an empty region
 * We define an empty region as follows:
 * region [i,j] is an empty region if for all k i<k<j, k is unpaired
 *
 * in our program we only need to know empty regions based on the original structure
 *
 * for every base i and j, such that i <= j
 * 1) if (i == j && i is not paired in G)
 * then region [i,i] is an empty region and
 * we put not_paired_all[i,j] = 1
 *
 * 2) else if j is not paired in G and we know that region[i,j-1] is an empty region
 * then region [i,j] is an empty region and we put
 * not_paired_all[i,j] = 1
 *
 * 3) otherwise region[i,j] is not an empty region and we put
 * not_paired_all[i,j] = 0
 *
 */

void detect_not_paired_all(h_str_features *fres, int *not_paired_all, int nb_nucleotides, int *index){
	int i, j;
	for(i = 0; i < nb_nucleotides; i++){
		for(j = i; j < nb_nucleotides; j++){
			int ij = index[i]+j-i;
			if (i == j){
				if (fres[i].pair < 0){
					not_paired_all[ij]=1;
//					if (debug){//_empty){
//						printf("region [%d,%d] is empty \n",i,j);
//					}
				}else{
					not_paired_all[ij] = 0;
//					if (debug){//_empty){
//						printf("region [%d,%d] is NOT empty \n",i,j);
//					}
				}
			}else if (not_paired_all[ij-1] == 1 && fres[j].pair < 0){
				not_paired_all[ij] = 1;
//				if (debug){//_empty){
//					printf("region [%d,%d] is empty \n",i,j);
//				}
			}
			else{
				not_paired_all[ij] = 0;
//				if (debug){//_empty){
//					printf("region [%d,%d] is NOT empty \n",i,j);
//				}
			}
		}
	}
}



/* Hosna:
 * The algorithm for finding b(i,l) and B(l,j) is as follows:
 *
 *  for every base l, if l is in an arc and l is not paired (i.e. f[l].arc != -1)
 *  and every i (0 < i < nb_nucleotides)
 *
 *  1) if (i <= arc(l)) then we are finding b(i,l):
 * 		1-a) temp = arc(arc(l))
 * 		1-b) if (temp == -1 || temp <= i) then we put
 * 			b(i,l) = b'(i,l) (i.e. border_bs[l][i] = arc(l))
 * 		1-c) otherwise
 * 				1-c-i) while (arc(temp) ! = -1 && arc(temp) > i)
 * 							temp = arc(temp)
 * 				1-c-ii) b(i,l) = temp and we put border_bs[l][i] = temp
 *
 *  2) if (i >= pair(arc(l))) then we are finding B(l,j):
 * 		1-a) temp = arc(arc(l))
 * 		1-b) if (temp == -1 || pair(temp) >= i) then we put
 * 			B(l,i) = B'(l,i) (i.e. border_bs[l][i] = pair(arc(l)))
 * 		1-c) otherwise
 * 				1-c-i) while (arc(temp) ! = -1 && pair(arc(temp)) < i)
 * 							temp = arc(temp)
 * 				1-c-ii) B(l,i) = pair(temp) and we put border_bs[l][i] = pair(temp)
 *
 *  3) if (i > arc(l) && i < pair(arc(l))) then border_bs[l][i] = -1
 */

void detect_border_bs(h_str_features *fres, int** border_bs, int nb_nucleotides){

	int l,i;
	for (l = 0; l < nb_nucleotides; l++){
		for (i = 0; i < nb_nucleotides; i++){
			int cover_l = fres[l].arc, pair_l=fres[l].pair; // Hosna March 8, 2012, using local varibales for optimality
			if (cover_l == -1 || pair_l >= 0){
				border_bs[l][i] = -2;
//				if (debug_border_bs){
//					printf("fres[%d].arc = %d, fres[%d].pair = %d \n",l, fres[l].arc,l,fres[l].pair);
//					printf("border_bs[%d][%d] = %d \n",l,i,border_bs[l][i]);
//				}

			}else{
				if (i <= cover_l){
					int temp = fres[cover_l].arc;
					if (temp == -1 || temp < i){ // Hosna: Jan 31, 2007: temp <= i changed to < to include i itself too
						border_bs[l][i] = cover_l;
					}else{
						while(fres[temp].arc != -1 && fres[temp].arc >= i){ // Hosna: Jan 31, 2007: fres[temp].arc > i changed to >= to include i itself too
							temp = fres[temp].arc;
						}
						border_bs[l][i] = temp;
					}
				}
				if (i >= fres[cover_l].pair){
					int temp = fres[cover_l].arc;
					if (temp == -1 || fres[temp].pair > i){ //Hosna: Jan 31, 2007: fres[temp].pair >= i changed to > to include i itself too
						border_bs[l][i] = fres[cover_l].pair;
					}else{
						while(fres[temp].arc != -1 && fres[fres[temp].arc].pair <= i){ //Hosna: Jan 31, 2007: fres[fres[temp].arc].pair < i changed to <= to include i itself too
							temp = fres[temp].arc;
						}
						border_bs[l][i] = fres[temp].pair;
					}
				}
				if (i > cover_l && i < fres[cover_l].pair){
					border_bs[l][i] = -1;
				}
//				if (debug_border_bs){
//					printf("fres[%d].arc = %d, fres[%d].pair = %d \n",l, fres[l].arc,l,fres[l].pair);
//					printf("border_bs[%d][%d] = %d \n",l,i,border_bs[l][i]);
//				}
			}
		}
	}
}

/* Hosna:
 * The algorithm for finding b'(i,l) and B'(l,j) is as follows:
 *
 *  for every base l, if l is in an arc and l is not paired (i.e. f[l].arc != -1)
 *  and every i (0 < i < nb_nucleotides)
 *
 * 	1) if (i < arc(l)) then b'(i,l) = arc(l)
 *  and we put the corresponding value in border_bps[l][i]
 *
 *  2) if (i > pair(arc(l))) then B'(l,i) = pair(arc(l))
 *  and we put the corresponding value in border_bps[l][i]
 *
 *  3) if (i >= arc(l) && i<= pair(arc(l))) then border[l][i] = -1
 *
 */

void detect_border_bps(h_str_features *fres, int** border_bps, int nb_nucleotides){
    // TODO ian should these be different if the i,l they are at is at linker (X) ?

	int l, i;
	for (l = 0; l < nb_nucleotides; l++){
		for (i = 0; i < nb_nucleotides ; i++){
			int cover_l = fres[l].arc, pair_l=fres[l].pair; // Hosna March 8, 2012, using local varibales for optimality
			if (cover_l  == -1 || pair_l >= 0){
				border_bps[l][i] = -2;//INF;//-2;
//				if (debug_border_bps){
//					printf("fres[%d].arc = %d, fres[%d].pair = %d \n",l, fres[l].arc,l,fres[l].pair);
//					printf("border_bps[%d][%d] = %d \n",l,i,border_bps[l][i]);
//				}
				//break;
			}else{
				if (i <= cover_l ){		//Hosna: Jan 31, 2007: < changed to <= to include i itself too
					border_bps[l][i] = cover_l ;
				}
				if (i >= fres[cover_l ].pair){ //Hosna: Jan 31, 2007: > changed to >= to include i itself too
					border_bps[l][i] = fres[fres[l].arc].pair;
				}
				if ( i > cover_l  && i < fres[cover_l].pair){ //Hosna: Jan 31, 2007: >= and <= changed to > and < to include i itself too
					border_bps[l][i] = -1;
				}
//				if (debug_border_bps){
//					printf("fres[%d].arc = %d, fres[%d].pair = %d \n",l, fres[l].arc,l,fres[l].pair);
//					printf("border_bps[%d][%d] = %d \n",l,i,border_bps[l][i]);
//				}
			}
		}
	}
}

void h_init (stack_ds *st)
// PRE:  None
// POST: Initialize the stack st
{
	st->top = STACK_EMPTY;//0;
}

void h_push (stack_ds *st, int el)
// PRE:  st is a valid stack
// POST: Push an element to the stack
{
	st->top = st->top +1;
    st->elem[st->top] = el;
}

int h_pop (stack_ds *st)
// PRE:  st is a valid stack, that is not empty
// POST: pop an element from the stack and return it
{
    if (st->top <= STACK_EMPTY)//0)
    {
        fprintf (stderr, "The given structure is not valid: more right parentheses than left parentheses\n");
        exit (1);
    }
    int result = st->elem[st->top];
    st->top = st->top -1 ;
    return result;
}

h_str_features *convert_str_features_to_h_str_features(str_features *f){
	h_str_features *fres;
	fres->arc = -1;
	fres->num_branches = f->num_branches;
	fres->pair = f->pair;
	fres->type = f->type;
	int i;
	for (i = 0; i < fres->num_branches; i++){
		fres->bri[i] = f->bri[i];
	}
	return fres;

}

str_features *convert_h_str_features_to_str_features(h_str_features *f){
	str_features *fres;
	fres->num_branches = f->num_branches;
	fres->pair = f->pair;
	fres->type = f->type;
	int i;
	for (i = 0; i < fres->num_branches; i++){
		fres->bri[i] = f->bri[i];
	}
	return fres;
}

void detect_h_structure_features (char *structure, h_str_features *f)
// PRE:  None
// POST: The variable f is filled with structure features, i.e. what type of elementary structure
//       this base is closing (such as stacked pair, hairpin loop etc.)
{
    int num_branches, i, j;
    int p_table[MAXSLEN];
    int arc_table[MAXSLEN];
    int bri[MAX_BRANCHES];
    int nb_nucleotides;
    nb_nucleotides = strlen(structure);
    detect_original_pairs_arcs (structure, p_table, arc_table);
    for (i=0; i < nb_nucleotides; i++)
    {
		// Hosna, March 8, 2012
		// use local variables instead of getting array values all the time
		int i_pair = p_table[i];
		f[i].pair = i_pair;//p_table[i];
        f[i].arc = arc_table[i];


        if (i_pair>i)//p_table[i] > i)
        {
            //f[i].pair = p_table[i]; //Hosna March 8, 2012, this statement seemed redundant! so removed
            f[i_pair].pair = i;//f[p_table[i]].pair = i;
            // check if it is stacked pair
			int i_pair_plus1 = p_table[i+1];
            if (i_pair_plus1 == i_pair-1 && i_pair_plus1 > i+1)//(p_table[i+1] == p_table[i]-1 && p_table[i+1] > i+1)
            {
                f[i].type = STACK;
                f[i_pair].type = STACK;//f[p_table[i]].type = STACK;
				/* if (i>=136){
						printf("END OF THE LOOP \n");
						printf("p_table[%d]=%d \n",i,p_table[i]);
						printf("f[%d].pair = %d, f[%d].type = %c, f[%d].arc = %d \n",i,f[i].pair,i,f[i].type,i,f[i].arc);

					} */
                continue;
            }
            // check if it is hairpin, internal loop or multi-loop
            num_branches = 0;
            for (j=i+1; j<i_pair; j++)//(j=i+1; j < p_table[i]; j++)
            {
                if (p_table[j] > j)
                {
                    bri[num_branches] = j;
                    num_branches++;
                    j = p_table[j];
                }
            }
            if (num_branches == 0)  // hairpin
            {
                f[i].type = HAIRP;
                f[i_pair].type = HAIRP; //f[p_table[i]].type = HAIRP;
            }
            else if (num_branches == 1) // internal loop
            {
                f[i].type = INTER;
                f[i_pair].type = INTER; //f[p_table[i]].type = INTER;
                f[i].num_branches = 1;
                f[i].bri[0] = bri[0];
            }
            else    // multi loop
            {
                f[i].type = MULTI;
                f[i_pair].type = MULTI; //f[p_table[i]].type = MULTI;
                f[i].num_branches = num_branches;
                for (j=0; j < num_branches; j++)
                    f[i].bri[j] = bri[j];
            }
        }

		/* if (i>=136){
							printf("END OF THE LOOP \n");
							printf("p_table[%d]=%d \n",i,p_table[i]);
							printf("f[%d].pair = %d, f[%d].type = %c, f[%d].arc = %d \n",i,f[i].pair,i,f[i].type,i,f[i].arc);

						} */
    }
    if (debug){
    	printf("h_str_features was successful! \n");
    }
}

/*
 * Hosna: January 10, 2008
 * The following two functions are modified versions of
 * the functions found in simfold/src/common/common.cpp
 * the modifications are to make them work for density-2 structures
 *
 */


double compute_h_sensitivity (char *ref_structure, char *pred_structure)
// returns 0 if undefined (denominator is 0)
{
    int ptable_ref[MAXSLEN];
    int ptable_pred[MAXSLEN];
    int distance;
    int len, i;
    double sens;
    int num_correct_bp;
    int num_true_bp;

    len = strlen(ref_structure);
    detect_original_PKed_pairs (ref_structure, ptable_ref);
    detect_original_PKed_pairs (pred_structure, ptable_pred);
    num_correct_bp = 0;
    num_true_bp = 0;
    for (i=0; i < len; i++)
    {
        if (ptable_ref[i] > -1)    // paired base
        {
            num_true_bp++;
            if (ptable_pred[i] == ptable_ref[i])
                num_correct_bp++;
        }
    }
    if (num_true_bp == 0)
        return -1.0;
    sens = num_correct_bp*1.0/num_true_bp;
    return sens;
}


double compute_h_ppv (char *ref_structure, char *pred_structure)
// returns 0 if undefined (denominator is 0)
{
    int ptable_ref[MAXSLEN];
    int ptable_pred[MAXSLEN];
    int distance;
    int len, i;
    double ppv;
    int num_correct_bp;
    int num_pred_bp;

    len = strlen(ref_structure);
    detect_original_PKed_pairs (ref_structure, ptable_ref);
    detect_original_PKed_pairs (pred_structure, ptable_pred);
    num_correct_bp = 0;
    num_pred_bp = 0;
    for (i=0; i < len; i++)
    {
        if (ptable_ref[i] > -1 && ptable_pred[i] == ptable_ref[i])    // paired base
            num_correct_bp++;
        if (ptable_pred[i] > -1)    // paired base
            num_pred_bp++;
    }
    if (num_pred_bp == 0)
        return -1.0;
    ppv = num_correct_bp*1.0/num_pred_bp;
    return ppv;
}

/*
 * Hosna Jan 10, 2008
 * The following function is used to tune the parameters
 * using Andronescu's GC program
 * This function is supposed to reset the pseudoknotted parameters
 */

void h_fill_data_structures_with_new_parameters (char *filename){
    FILE *file;
    char buffer[100];
    double param;
    int line = 0;

    //printf ("FILENAME: %s\n", filename);
	if ((file = fopen (filename, "r")) == NULL)
	{
	    giveup ("Cannot open file", filename);
	}

	// PS_penalty: exterior pseudoloop initiation penalty (originally 9.6 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
    param *= 100;
    PS_penalty = (int) param;

	//PSM_penalty: penalty for introducing pseudoknot inside a multiloop (originally 15 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    PSM_penalty = (int) param;

	//PSP_penalty: penalty for introducing pseudoknot inside a pseudoloop (originally 15 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    PSP_penalty = (int) param;

	//PB_penalty: band penalty (originally 0.2 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    PB_penalty = (int) param;

	//PUP_penalty: penalty for an un-paired base in a pseudoloop or a band (originally 0.1 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    PUP_penalty = (int) param;

	//PPS_penalty: penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band(originally 0.1 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    PPS_penalty = (int) param;


	//a_penalty: penalty for introducing a multiloop (originally 3.4 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    a_penalty = (int) param;

	//b_penalty: penalty for base pair in a multiloop (originally 0.4 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    b_penalty = (int) param;


	//c_penalty: penalty for un-paired base in a multi-loop (originally 0)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
    param *= 100;
    c_penalty = (int) param;



	// e_stP = 0.83 * e_s
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
    e_stP_penalty = param;


	// e_intP = 0.83 * e_int
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
    e_intP_penalty = param;


	//ap_penalty: penalty for introducing a multiloop that spans a band (originally 3.4 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    ap_penalty = (int) param;

	//bp_penalty: base pair penalty for a multiloop that spans a band (originally 0.4 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    bp_penalty = (int) param;


	//cp_penalty: penalty for unpaired base in a multiloop that spans a band (originally 0)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    cp_penalty = (int) param;


	fclose (file);
	//printf ("****** we must have 14 lines by now: LINES = %d\n", line);


}


void h_fill_data_structures_with_new_parameters (double *param)
{

	if (param == NULL){
		giveup ("Incorrect parameter length", "h_fill_data_structure_with_new_parameters");
	}

	// PS_penalty: exterior pseudoloop initiation penalty (originally 9.6 Kcal/mol)
    param[0] *= 100;
    PS_penalty = (int) param[0];

	//PSM_penalty: penalty for introducing pseudoknot inside a multiloop (originally 15 Kcal/mol)
	param[1] *= 100;
    PSM_penalty = (int) param[1];

	//PSP_penalty: penalty for introducing pseudoknot inside a pseudoloop (originally 15 Kcal/mol)
	param[2] *= 100;
    PSP_penalty = (int) param[2];

	//PB_penalty: band penalty (originally 0.2 Kcal/mol)
	param[3] *= 100;
    PB_penalty = (int) param[3];

	//PUP_penalty: penalty for an un-paired base in a pseudoloop or a band (originally 0.1 Kcal/mol)
	param[4] *= 100;
    PUP_penalty = (int) param[4];

	//PPS_penalty: penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band(originally 0.1 Kcal/mol)
    param[5] *= 100;
    PPS_penalty = (int) param[5];


	//a_penalty: penalty for introducing a multiloop (originally 3.4 Kcal/mol)
	param[6] *= 100;
    a_penalty = (int) param[6];

	//b_penalty: penalty for base pair in a multiloop (originally 0.4 Kcal/mol)
    param[7] *= 100;
    b_penalty = (int) param[7];


	//c_penalty: penalty for un-paired base in a multi-loop (originally 0)
    param[8] *= 100;
    c_penalty = (int) param[8];



	// e_stP = 0.83 * e_s
    e_stP_penalty = param[9];


	// e_intP = 0.83 * e_int
    e_intP_penalty = param[10];


	//ap_penalty: penalty for introducing a multiloop that spans a band (originally 3.4 Kcal/mol)
    param[11] *= 100;
    ap_penalty = (int) param[11];

	//bp_penalty: base pair penalty for a multiloop that spans a band (originally 0.4 Kcal/mol)
    param[12] *= 100;
    bp_penalty = (int) param[12];


	//cp_penalty: penalty for unpaired base in a multiloop that spans a band (originally 0)
    param[13] *= 100;
    cp_penalty = (int) param[13];

}

int h_create_string_params(){
	int index = create_string_params();

	sprintf (string_params[index++], "PS_penalty");
	sprintf (string_params[index++], "PSM_penalty");
	sprintf (string_params[index++], "PSP_penalty");
	sprintf (string_params[index++], "PB_penalty");
	sprintf (string_params[index++], "PUP_penalty");
	sprintf (string_params[index++], "PPS_penalty");


	sprintf (string_params[index++], "a_penalty");
	sprintf (string_params[index++], "b_penalty");
	sprintf (string_params[index++], "c_penalty");


	sprintf (string_params[index++], "e_stP_penalty");
	sprintf (string_params[index++], "e_intP_penalty");


	sprintf (string_params[index++], "ap_penalty");
	sprintf (string_params[index++], "bp_penalty");
	sprintf (string_params[index++], "cp_penalty");
	return index;
}

double hfold(char *sequence, char *restricted, char *structure){
	W_final *min_fold = new W_final (sequence, restricted);
	if (min_fold == NULL) giveup ("Cannot allocate memory", "HFold");

	double energy = min_fold->hfold();
    min_fold->return_structure (structure);

    delete min_fold;
    return energy;
}

//AP. HFold_Emodel is used by HFold and HFold_interacting. The only difference is the passed in energy model vector.
/*
double hfold_emodel(char *sequence, char *restricted, char *structure, std::vector<energy_model> *energy_models){
	W_final *min_fold = new W_final (sequence, restricted, energy_models);
	if (min_fold == NULL) giveup ("Cannot allocate memory", "HFold");
	double energy = min_fold->hfold_emodel();
    min_fold->return_structure (structure);
    return energy;
}
*/

//---------------------------------------this function is suppose to be the same as the one in Hfold_iterative, if any changes are made, please change that one too--------------------
//kevin 19 july (verified by Mahyar 19 july 2017)
//count+1 when open bracket, count-1 when close bracket
//whgen count is 0, we have a substructure
//assume valid structure
void find_disjoint_substructure(char* structure, std::vector< std::pair<int,int> > &pair_vector){
	int length = strlen(structure);
	int count = 0;
	int first_time = 1; //flag for first time getting a open bracket for substructure
	int i = 0;
	int j = 0;
	for(int k=0; k<length;k++){
		if(structure[k] == '(' || structure[k] == '['){
			if(first_time && count == 0){
				first_time = 0;
				i = k;
			}
			count += 1;

		}else if(structure[k] == ')' || structure[k] == ']'){
			count -= 1;
			j = k;
			if(count == 0){
				std::pair <int,int> ij_pair (i,j);
				pair_vector.push_back(ij_pair);
				first_time = 1;
			}
		}
	}
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_iterative, if any changes are made, please change that one too--------------------
//31 Aug 2017 kevin and Mahyar
//input[i] is _ or . and it did not turn into a . in the output structure, then it is not empty
int is_empty_structure(char* input_structure, char* output_structure){
	for(int i=0; i<strlen(input_structure);i++){
		if(input_structure[i] != output_structure[i]){
			if((input_structure[i] == '_' || input_structure[i] == '.' ) && output_structure[i] != '.'){
				return 0;
			}
		}
	}
	return 1;
}




//---------------------------------------this function is suppose to be the same as the one in Hfold_iterative, if any changes are made, please change that one too--------------------
//kevin 18 July
int paired_structure(int i, int j, int *pair_index, int length){
	if(i >= 0 && j < length && (pair_index[i] == j) && (pair_index[j] == i) ){
		return 1;
	}
	return 0;
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_iterative, if any changes are made, please change that one too--------------------
//kevin 18 July
void obtainRelaxedStems(char* G1, char* G2, char* Gresult){
	int length = strlen(G1);
	int G1_pair[length];
	int G2_pair[length];

	//Gresult <- G1
	strncpy(Gresult,G1,length);

	detect_original_pairs(G1, G1_pair);
	detect_original_pairs(G2, G2_pair);

	//for(int d=0;d<length;d++){
	//	printf("%c %d\n",G2[d],G2_pair[d]);
	//}

	int i = 0;
	int j = 0;

	for(int k=0;k<length;k++){
		if(G2_pair[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (G1[i] != G2[i]) && (G1[j] != G2[j]) ){//if ij not in G1
					//include bulges of size 1
					if(paired_structure(i-1,j+1,G1_pair,length) || paired_structure(i+1,j-1,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					//include loops of size 1x1
					}else if( paired_structure(i-2,j+1,G1_pair,length) || paired_structure(i-1,j+2,G1_pair,length) || \
							paired_structure(i+1,j-2,G1_pair,length) || paired_structure(i+2,j-1,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					//include loops of size 1x2 or 2x1
					}else if( paired_structure(i-2,j+2,G1_pair,length) || paired_structure(i+2,j-2,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					}else if( paired_structure(i-3,j+2,G1_pair,length) || paired_structure(i-2,j+3,G1_pair,length) || \
							paired_structure(i+2,j-3,G1_pair,length) || paired_structure(i+3,j-2,G1_pair,length) ){

						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					}
				}
			}
		}
	}
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_iterative, if any changes are made, please change that one too--------------------
//kevin 30 Aug 2017
//check if the computed structure matches the restricted structure
int is_invalid_restriction(char* restricted_structure, char* current_structure){
	for (int i=0; i < strlen(restricted_structure); i++){
        if ((restricted_structure[i] == '(' || restricted_structure[i] == ')' || restricted_structure[i] == '.') &&
            (restricted_structure[i] != current_structure[i])){
				return 1;
        }
    }
	return 0;
}


//kevin 18 July
void simfold_emodel(char *sequence, char *restricted, char *structure, std::vector<energy_model> *energy_models){
	W_final *simfold = new W_final (sequence, restricted, energy_models);
	if (simfold == NULL) giveup ("Cannot allocate memory", "method3 Simfold");
	double min_energy = simfold->call_simfold_emodel();
	simfold->return_structure (structure);

    delete simfold;
}

//kevin 24 July
double method1_emodel(char *sequence, char *restricted, char *structure, std::vector<energy_model> *energy_models){
	W_final *hfold_min_fold = new W_final (sequence, restricted, energy_models);
	if (hfold_min_fold == NULL) giveup ("Cannot allocate memory", "HFold");
	double energy = 0;
	energy = hfold_min_fold->hfold_emodel();
	hfold_min_fold->return_structure (structure);
	
    delete hfold_min_fold;
	return energy;
}


//kevin 24 July
double method2_emodel(char *sequence, char *restricted, char *structure, std::vector<energy_model> *energy_models){
	double energy = 0;
	W_final *hfold_pk_min_fold = new W_final (sequence, restricted, energy_models);
	if (hfold_pk_min_fold == NULL) giveup ("Cannot allocate memory", "HFold");
	energy = hfold_pk_min_fold->hfold_pkonly();
	hfold_pk_min_fold->return_structure (structure);

	if(is_empty_structure(restricted,structure)){
    	delete hfold_pk_min_fold;
		return energy;
	}else{
		char G_prime[strlen(structure)];
		remove_structure_intersection(structure,restricted, G_prime);

		energy = method1_emodel(sequence, G_prime, structure, energy_models);

    	delete hfold_pk_min_fold;
		return energy;
	}

}

//kevin 18 July
double method3_emodel(char *sequence, char *restricted, char *structure, std::vector<energy_model> *energy_models){
	double energy = 0;
	int length = strlen(sequence);
	char simfold_structure[length];

	simfold_emodel(sequence,restricted, simfold_structure, energy_models);
	//^ G' simfold_structure <- SimFold(S sequence, G restricted)
	char G_updated[length+1];
	G_updated[length] = '\0';
	obtainRelaxedStems(restricted ,simfold_structure, G_updated);
	//^Gupdated G_updated<- ObtainRelaxedStems(G restricted,G' simfold_structure)
	energy = method2_emodel(sequence, G_updated, structure, energy_models);

	return energy;
}

//kevin 18 July
double method4_emodel(char *sequence, char *restricted, char *structure, std::vector<energy_model> *energy_models){
	int KEVIN_DEBUG = 0;
	double energy = 0;
	int length = strlen(sequence);
	char G_updated[length+1];
	int k = 1;
	//^k <- 1
	strcpy(G_updated, restricted);
	//^Gupdated <- G
	std::vector< std::pair<int,int> > disjoint_substructure_index; //contain pair(i,j) in a vector, ij are index of begin and end of substructure
	find_disjoint_substructure(restricted, disjoint_substructure_index);
	//^get disjoint substructure
	int i = 0;
	int j = 0;
	for(auto current_substructure_index : disjoint_substructure_index){
		i = current_substructure_index.first;
		j = current_substructure_index.second;
		char subsequence[length+1];
		char substructure[length+1];
		char simfold_structure[length+1];

		strncpy(subsequence, sequence+i,j-i+1);
		subsequence[j-i+1] = '\0';
		//^Sk
		strncpy(substructure, restricted+i,j-i+1);
		substructure[j-i+1] = '\0';
		//^Gk

        // Now that we have split into substructure, linker has moved
        // so keep older linker location for restoring after
        int temp_linker_pos = linker_pos;

        // if whole subsequence is left of linker_pos,
        if((i < linker_pos) && (j < linker_pos)){ //left of linker
            // set linker_pos to INF. This way subsequence is always considered to left of linker
            linker_pos = INF;
        }
        // if whole subsequence is right of linker_pos
        else if((i > linker_pos+linker_length-1) && (j > linker_pos+linker_length-1)){ //right of linker
            // set linker_pos to linker_length.
            // can't be 0 because ((i > linker_pos+linker_length-1) && (j > linker_pos+linker_length-1)) means any i or j below linker_length incorrectly returns false
            // This way subsequence is always considered to right of linker
            linker_pos = -linker_length;
        }
        // if subsequence contains linker
        else if((i < linker_pos) && (j > linker_pos+linker_length-1)){ //cross linker
            // linker_pos is position of beginning of linker, so just move it over by the beginning of our subsequence
            linker_pos = linker_pos-i;
        } else {
            fprintf(stderr,"ERROR in method 4 setting new linker_pos\n");
        }
		simfold_emodel(subsequence, substructure, simfold_structure, energy_models);
		//^ SimFold(Sk,Gk,Gk',energy_models)
		char Gp_k_updated[length];
		obtainRelaxedStems(substructure, simfold_structure, Gp_k_updated);
		//^obtainRelaxedStems(Gk,Gk',G'kupdated)
		linker_pos = temp_linker_pos;
		int m = 0; //index for going through Gp_k_updated
        for(int k =i;k<j;k++){
			if(G_updated[k] != Gp_k_updated[m]){
				G_updated[k] = Gp_k_updated[m];
			}
			m++;
		}
		//^Gupdated <- Gupdated U G'kupdated
	}
	energy = method2_emodel(sequence, G_updated, structure, energy_models);

	return energy;
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_iterative, if any changes are made, please change that one too--------------------
// Aug 31, 2017 kevin and Mahyar
//does G_p = G1-G
void remove_structure_intersection(char* G1, char* G, char* G_p){
	strcpy(G_p,G1);
	for(int i=0; i< strlen(G1); i++){  
		if (G1[i] != G[i]){
			continue;
		}else{
			G_p[i] = '.';
		}
	}
}


//kevin 18 July
//july 24: changed hfold, hfold_pkonly to a method; changed replaced final_structure with method1-4_structure
double hfold_interacting_emodel(char *sequence, char *restricted, char *structure, std::vector<energy_model> *energy_models, int &method_used){

	double energy = 0;
	double min_energy = INF;
	char method1_structure[strlen(sequence)+1];
	char method2_structure[strlen(sequence)+1];
	char method3_structure[strlen(sequence)+1];
	char method4_structure[strlen(sequence)+1];
printf("start method1\n");
	min_energy = method1_emodel(sequence,restricted,method1_structure,energy_models);
	method_used = 1;
	strcpy(structure,method1_structure);
	//printf("method1 energy: %lf\n",min_energy);

printf("start method2\n");
	energy = method2_emodel(sequence,restricted,method2_structure,energy_models);
	if(energy < min_energy){
        method_used = 2;
		min_energy = energy;
		strcpy(structure,method2_structure);
	}

printf("start method3\n");
	energy = method3_emodel(sequence,restricted,method3_structure,energy_models);
	if(energy < min_energy){
        method_used = 3;
		min_energy = energy;
		strcpy(structure,method3_structure);
    }
printf("start method4\n");
	energy = method4_emodel(sequence,restricted,method4_structure,energy_models);
	if(energy < min_energy){
        method_used = 4;
		min_energy = energy;
		strcpy(structure,method4_structure);
	}

	return min_energy;
}


// Hosna, May 3rd, 2012
// added function for the pkonly version
double hfold_pkonly(char *sequence, char *restricted, char *structure){
	W_final *min_fold = new W_final (sequence, restricted);
	if (min_fold == NULL) giveup ("Cannot allocate memory", "HFoldPKonly");

	double energy = min_fold->hfold_pkonly();
    min_fold->return_structure (structure);

    delete min_fold;
    return energy;
}

//AP. HFold_Emodel is used by HFold_pkonly and HFold_interacting_pkonly. The only difference is the passed in energy model vector.
double hfold_pkonly_emodel(char *sequence, char *restricted, char *structure, std::vector<energy_model> *energy_models){
	W_final *min_fold = new W_final (sequence, restricted, energy_models);
	if (min_fold == NULL) giveup ("Cannot allocate memory", "HFoldPKonly");

	double energy = min_fold->hfold_pkonly_emodel();
    min_fold->return_structure (structure);

    delete min_fold;
    return energy;
}



//AP. This function is used to print out the data inside an energy model. Used for debugging only.
void print_emodel(energy_model *model) {
	FILE *ioFile;
	char filePath[MAXSLEN];

	/*printf("similarity_rule\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "similarity_rule.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXNUMPARAMS; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < 5000; col2++) {
			fprintf(ioFile, "%c ", model->similarity_rule[col1][col2]);
			//printf("%c ", model->similarity_rule[col1][col2]);
			fflush(stdout);
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);	*/

	printf("\nstack\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "stack.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "%d", model->stack[col1][col2][col3][col4]);
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\ntstackh\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "tstackh.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "%d", model->tstackh[col1][col2][col3][col4]);
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\ntstacki\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "tstacki.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "%d", model->tstacki[col1][col2][col3][col4]);
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nint11\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "int11.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "{");
					for(int col5 = 0; col5 < NUCL; col5++) {
						fprintf(ioFile, "{");
						for(int col6 = 0; col6 < NUCL; col6++) {
							fprintf(ioFile, "%d", model->int11[col1][col2][col3][col4][col5][col6]);
						}
						fprintf(ioFile, "}");
					}
					fprintf(ioFile, "}");
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nint21\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "int21.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "{");
					for(int col5 = 0; col5 < NUCL; col5++) {
						fprintf(ioFile, "{");
						for(int col6 = 0; col6 < NUCL; col6++) {
							fprintf(ioFile, "{");
							for(int col7 = 0; col7 < NUCL; col7++) {
								fprintf(ioFile, "%d", model->int21[col1][col2][col3][col4][col5][col6][col7]);
							}
							fprintf(ioFile, "}");
						}
						fprintf(ioFile, "}");
					}
					fprintf(ioFile, "}");
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nint22\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "int22.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "{");
					for(int col5 = 0; col5 < NUCL; col5++) {
						fprintf(ioFile, "{");
						for(int col6 = 0; col6 < NUCL; col6++) {
							fprintf(ioFile, "{");
							for(int col7 = 0; col7 < NUCL; col7++) {
								fprintf(ioFile, "{");
								for(int col8 = 0; col8 < NUCL; col8++) {
									fprintf(ioFile, "%d", model->int22[col1][col2][col3][col4][col5][col6][col7][col8]);
								}
								fprintf(ioFile, "}");
							}
							fprintf(ioFile, "}");
						}
						fprintf(ioFile, "}");
					}
					fprintf(ioFile, "}");
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\ndangle_top\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "dangle_top.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "%d", model->dangle_top[col1][col2][col3]);
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\ndangle_bot\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "dangle_bot.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "%d", model->dangle_bot[col1][col2][col3]);
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\ninternal_penalty_by_size\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "internal_penalty_by_size.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXLOOP+1; col1++) {
		fprintf(ioFile, "%d", model->internal_penalty_by_size[col1]);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nbulge_penalty_by_size\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "bulge_penalty_by_size.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXLOOP+1; col1++) {
		fprintf(ioFile, "%d", model->bulge_penalty_by_size[col1]);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nhairpin_penalty_by_size\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "hairpin_penalty_by_size.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXLOOP+1; col1++) {
		fprintf(ioFile, "%d", model->hairpin_penalty_by_size[col1]);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\ntriloop\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "triloop.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXTRILOOPNO; col1++) {
		fprintf(ioFile, "{%s%d}", model->triloop[col1].seq, model->triloop[col1].energy);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\ntloop\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "tloop.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXTRILOOPNO; col1++) {
		fprintf(ioFile, "{%s%d}", model->tloop[col1].seq, model->tloop[col1].energy);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nnb_triloops\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "nb_triloops.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "%d", model->nb_triloops);
	fclose(ioFile);

	printf("\nnb_tloops\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "nb_tloops.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "%d", model->nb_tloops);
	fclose(ioFile);

	printf("\nspecial_hl\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "special_hl.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAX_SPECIAL_LOOP_NO; col1++) {
		fprintf(ioFile, "{%s%d}", model->special_hl[col1].seq, model->special_hl[col1].energy);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nnb_special_hl\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "nb_special_hl.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "%d", model->nb_special_hl);
	fclose(ioFile);

	printf("\nint11_experimental_addition\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "int11_experimental_addition.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "{");
					for(int col5 = 0; col5 < NUCL; col5++) {
						fprintf(ioFile, "{");
						for(int col6 = 0; col6 < NUCL; col6++) {
							fprintf(ioFile, "%d", model->int11_experimental_addition[col1][col2][col3][col4][col5][col6]);
						}
						fprintf(ioFile, "}");
					}
					fprintf(ioFile, "}");
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "z}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nint21_experimental_addition\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "int21_experimental_addition.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "{");
					for(int col5 = 0; col5 < NUCL; col5++) {
						fprintf(ioFile, "{");
						for(int col6 = 0; col6 < NUCL; col6++) {
							fprintf(ioFile, "{");
							for(int col7 = 0; col7 < NUCL; col7++) {
								fprintf(ioFile, "%d", model->int21_experimental_addition[col1][col2][col3][col4][col5][col6][col7]);
							}
							fprintf(ioFile, "}");
						}
						fprintf(ioFile, "}");
					}
					fprintf(ioFile, "}");
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nint22_experimental_addition\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "int22_experimental_addition.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "{");
					for(int col5 = 0; col5 < NUCL; col5++) {
						fprintf(ioFile, "{");
						for(int col6 = 0; col6 < NUCL; col6++) {
							fprintf(ioFile, "{");
							for(int col7 = 0; col7 < NUCL; col7++) {
								fprintf(ioFile, "{");
								for(int col8 = 0; col8 < NUCL; col8++) {
									fprintf(ioFile, "%d", model->int22_experimental_addition[col1][col2][col3][col4][col5][col6][col7][col8]);
								}
								fprintf(ioFile, "}");
							}
							fprintf(ioFile, "}");
						}
						fprintf(ioFile, "}");
					}
					fprintf(ioFile, "}");
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nbulge1\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "bulge1.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "{");
					for(int col5 = 0; col5 < NUCL; col5++) {
						fprintf(ioFile, "%d", model->bulge1[col1][col2][col3][col4][col5]);
					}
					fprintf(ioFile, "}");
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\ninternal_asymmetry\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "internal_asymmetry.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXLOOP+1; col1++) {
		fprintf(ioFile, "%d", model->internal_asymmetry[col1]);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_stack\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_stack.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "%d", model->enthalpy_stack[col1][col2][col3][col4]);
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_tstackh\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_tstackh.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "%d", model->enthalpy_tstackh[col1][col2][col3][col4]);
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_tstacki\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_tstacki.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "%d", model->enthalpy_tstacki[col1][col2][col3][col4]);
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_int11\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_int11.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "{");
					for(int col5 = 0; col5 < NUCL; col5++) {
						fprintf(ioFile, "{");
						for(int col6 = 0; col6 < NUCL; col6++) {
							fprintf(ioFile, "%d", model->enthalpy_int11[col1][col2][col3][col4][col5][col6]);
						}
						fprintf(ioFile, "}");
					}
					fprintf(ioFile, "}");
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_int21\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_int21.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "{");
					for(int col5 = 0; col5 < NUCL; col5++) {
						fprintf(ioFile, "{");
						for(int col6 = 0; col6 < NUCL; col6++) {
							fprintf(ioFile, "{");
							for(int col7 = 0; col7 < NUCL; col7++) {
								fprintf(ioFile, "%d", model->enthalpy_int21[col1][col2][col3][col4][col5][col6][col7]);
							}
							fprintf(ioFile, "}");
						}
						fprintf(ioFile, "}");
					}
					fprintf(ioFile, "}");
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_int22\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_int22.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "{");
				for(int col4 = 0; col4 < NUCL; col4++) {
					fprintf(ioFile, "{");
					for(int col5 = 0; col5 < NUCL; col5++) {
						fprintf(ioFile, "{");
						for(int col6 = 0; col6 < NUCL; col6++) {
							fprintf(ioFile, "{");
							for(int col7 = 0; col7 < NUCL; col7++) {
								fprintf(ioFile, "{");
								for(int col8 = 0; col8 < NUCL; col8++) {
									fprintf(ioFile, "%d", model->enthalpy_int22[col1][col2][col3][col4][col5][col6][col7][col8]);
								}
								fprintf(ioFile, "}");
							}
							fprintf(ioFile, "}");
						}
						fprintf(ioFile, "}");
					}
					fprintf(ioFile, "}");
				}
				fprintf(ioFile, "}");
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_dangle_top\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_dangle_top.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "%d", model->enthalpy_dangle_top[col1][col2][col3]);
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_dangle_bot\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_dangle_bot.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < NUCL; col1++) {
		fprintf(ioFile, "{");
		for(int col2 = 0; col2 < NUCL; col2++) {
			fprintf(ioFile, "{");
			for(int col3 = 0; col3 < NUCL; col3++) {
				fprintf(ioFile, "%d", model->enthalpy_dangle_bot[col1][col2][col3]);
			}
			fprintf(ioFile, "}");
		}
		fprintf(ioFile, "}");
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_internal_penalty_by_size\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_internal_penalty_by_size.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXLOOP+1; col1++) {
		fprintf(ioFile, "%d", model->enthalpy_internal_penalty_by_size[col1]);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_bulge_penalty_by_size\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_bulge_penalty_by_size.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXLOOP+1; col1++) {
		fprintf(ioFile, "%d", model->enthalpy_bulge_penalty_by_size[col1]);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_hairpin_penalty_by_size\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_hairpin_penalty_by_size.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXLOOP+1; col1++) {
		fprintf(ioFile, "%d", model->enthalpy_hairpin_penalty_by_size[col1]);
	}

	printf("\nenthalpy_triloop\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_triloop.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXTRILOOPNO; col1++) {
		fprintf(ioFile, "{%s%d}", model->enthalpy_triloop[col1].seq, model->enthalpy_triloop[col1].energy);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_tloop\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_tloop.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "{");
	for(int col1 = 0; col1 < MAXTRILOOPNO; col1++) {
		fprintf(ioFile, "{%s%d}", model->enthalpy_tloop[col1].seq, model->enthalpy_tloop[col1].energy);
	}
	fprintf(ioFile, "}");
	fclose(ioFile);

	printf("\nenthalpy_nb_triloops\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_nb_triloops.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "%d", model->enthalpy_nb_triloops);
	fclose(ioFile);

	printf("\nenthalpy_nb_tloops\n");
	strcpy(filePath, "./temp/");
	strcat(filePath, "enthalpy_nb_tloops.txt");
	ioFile = fopen(filePath, "w");
	fprintf(ioFile, "%d", model->enthalpy_nb_tloops);
	fclose(ioFile);
}
