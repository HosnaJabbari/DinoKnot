/***************************************************************************
                          s_energy_matrix.cpp  -  description
                             -------------------
    begin                : Fri Apr 12 2002
    copyright            : (C) 2002 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

 // This is the V matrix

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
// Hosna, March 5, 2012
// malloc.h is not needed in my mac as stdlib.h does the same
//#include <malloc.h>

#include "constants.h"
#include "structs.h"
#include "externs.h"
#include "common.h"
#include "simfold.h"
#include "s_energy_matrix.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"

#include "../src/pseudo_loop.h"

s_energy_matrix::s_energy_matrix (int *seq, int length)
// The constructor
{
    this->H = NULL;
    this->S = NULL;
    this->VBI = NULL;
    this->VM = NULL;

    sequence = seq;     // just refer it from where it is in memory
    seqlen = length;

    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    index = new int [length];
    int total_length = (length *(length+1))/2;
    index[0] = 0;
    for (int i=1; i < length; i++)
        index[i] = index[i-1]+length-i+1;

    // this array holds V(i,j), and what (i,j) encloses: hairpin loop, stack pair, internal loop or multi-loop
    nodes = new free_energy_node [total_length];
    if (nodes == NULL) giveup ("Cannot allocate memory", "s_energy_matrix");
}

//AP. Needed a custom constructor to add the energy moodel vector.
s_energy_matrix::s_energy_matrix (int *seq, int length, std::vector<energy_model> *energy_models)
// The constructor
{
	this->energy_models = energy_models;
    this->H = NULL;
    this->S = NULL;
    this->VBI = NULL;
    this->VM = NULL;

    sequence = seq;     // just refer it from where it is in memory
    seqlen = length;

    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    index = new int [length];
    int total_length = (length *(length+1))/2;
    index[0] = 0;
    for (int i=1; i < length; i++)
        index[i] = index[i-1]+length-i+1;

    // this array holds V(i,j), and what (i,j) encloses: hairpin loop, stack pair, internal loop or multi-loop
    nodes = new free_energy_node [total_length];
    if (nodes == NULL) giveup ("Cannot allocate memory", "s_energy_matrix");
}

s_energy_matrix::~s_energy_matrix ()
// The destructor
{
    delete [] index;
    delete [] nodes;
}



void s_energy_matrix::compute_energy (int i, int j)
// compute the V(i,j) value
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

    if (can_pair (sequence[i], sequence[j]))    // if i and j can pair
    {
        // compute free energy of hairpin loop, stack pair, internal loop and multi-loop
        min_en[0] = H->compute_energy (i, j);
        if (i<=j-TURN-1)
        {
            min_en[1] = S->compute_energy (i, j);

            // TODO: uncomment
            if (!ignore_internal)
                min_en[2] = VBI->compute_energy (i, j);
            if (!ignore_multi)
                min_en[3] = VM->compute_energy (i, j);
        }
    }

    // see which of them is the minimum
    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = STACK; break;
        case  2: type = INTER; break;
        case  3: type = MULTI; break;
        default: type = NONE;
    }

    if (min_rank > -1)
    {
        if (debug)
            printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
    }

    if (min < INF/2)
    {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}


void s_energy_matrix::compute_energy_restricted (int i, int j, str_features *fres)
// compute the V(i,j) value, if the structure must be restricted
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

	// Hosna, March 26, 2012
	// if the restricted base pairs are non-canonical then checking for can_pair only will cause missing those base pairs
	if (can_pair (sequence[i], sequence[j]) || (fres[i].pair == j && fres[j].pair ==i))
    {
        if (fres[i].pair == i+1)
            min_en[0] = 0;
        else
        {
            if (!exists_restricted (i, j, fres))
                min_en[0] = H->compute_energy_restricted (i, j, fres); // there was a stupid bug here, I was calling H->compute_energy instead of the restricted version. Fixed on June 30, 2007.

            min_en[1] = S->compute_energy_restricted (i, j,fres);//S->compute_energy (i, j); Hosna, March 26, 2012
            min_en[2] = VBI->compute_energy_restricted (i, j, fres);
            min_en[3] = VM->compute_energy_restricted (i, j, fres);
        }
    }

    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = STACK; break;
        case  2: type = INTER; break;
        case  3: type = MULTI; break;
        default: type = NONE;
    }

    if (min_rank > -1 && debug) {
    	printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
   	}

    if (min < INF/2) {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}

//AP
void s_energy_matrix::compute_energy_restricted_emodel (int i, int j, str_features *fres, int is_weakly_closed_flag) 
// compute the V(i,j) value, if the structure must be restricted
{
    
    PARAMTYPE min, min_en[4];
    int k, min_rank = -1;
    char type;

    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;
   
	// Hosna, March 26, 2012
	// if the restricted base pairs are non-canonical then checking for can_pair only will cause missing those base pairs
	if (can_pair (sequence[i], sequence[j]) || (fres[i].pair == j && fres[j].pair ==i)) {
        if (fres[i].pair == i+1)
            min_en[0] = 0;
        else
        {
			if (!exists_restricted (i, j, fres)) {
				for (auto &energy_model : *energy_models) {
                    //kevin question: this is one model, we check if type are different and add it here?
		            energy_model.energy_value = H->compute_energy_restricted_emodel (i, j, fres, &energy_model); // there was a stupid bug here, I was calling H->compute_energy instead of the restricted version. Fixed on June 30, 2007.
				}
				min_en[0] = emodel_energy_function (i, j, energy_models);

                // Ian Wark and Kevin July 20 2017
                // Hybrid molecule penalty
                // This needs to be out here as we cannot tell if the energy models are different in H->compute_energy_restricted_emodel
                // requires there to be exactly 2 energy models
				if( is_cross_model(i,j) ) {    // If cross model 
                   min_en[0] += START_HYBRID_PENALTY;  
                   //printf("added penalty\n");                                              // add a hybrid molecule penalty
				}
			}
            //16 Aug 2017 kevin
            //added if is_weakly_closed_flag to make sure there are no base pair going out of region ij, aka is_weakly_closed
            //have it as a flag instead of doing the check here because is_weakly_closed(int i, int j) is a method from the class pseudo_loop which is out of this class
            //so we have to do the check before we call this function which is W_final, a class containing a pseudo_loop object
            if(is_weakly_closed_flag){ 
                for (auto &energy_model : *energy_models) {
                    //S->compute_energy (i, j); Hosna, March 26, 2012
                    energy_model.energy_value = S->compute_energy_restricted_emodel (i, j,fres, &energy_model);
                }
                min_en[1] = emodel_energy_function (i, j, energy_models);

                for (auto &energy_model : *energy_models) {
                    energy_model.energy_value = VBI->compute_energy_restricted_emodel (i, j, fres, &energy_model);
                }
                min_en[2] = emodel_energy_function (i, j, energy_models);
                

                for (auto &energy_model : *energy_models) {
                    energy_model.energy_value = VM->compute_energy_restricted_emodel (i, j, fres, &energy_model);
                }
                
                min_en[3] = emodel_energy_function (i, j, energy_models);

            }
        }
    }

    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            
            min = min_en[k];
            min_rank = k;
        }
    }

    


    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = STACK; break;
        case  2: type = INTER; break;
        case  3: type = MULTI; break;
        default: type = NONE;
    }

//printf ("V(%d,%d) %c energy %d %d %d %d, min=%d\n", i, j, type, min_en[0],min_en[1],min_en[2],min_en[3],min);

    if (min < INF/2) {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }

}

//AP
void s_energy_matrix::compute_energy_restricted_pkonly_emodel (int i, int j, str_features *fres, int is_weakly_closed_flag)
// compute the V(i,j) value, if the structure must be restricted
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

	if (fres[i].pair == j && fres[j].pair ==i) {
        if (fres[i].pair == i+1) {
            min_en[0] = 0;
        } else {
			if (!exists_restricted (i, j, fres)) {
				for (auto &energy_model : *energy_models) {
	            	energy_model.energy_value = H->compute_energy_restricted_emodel (i, j, fres, &energy_model);
				}
				min_en[0] = emodel_energy_function (i, j, energy_models);

				// Ian Wark and Kevin July 20 2017
                // Hybrid molecule penalty
                // This needs to be out here as we cannot tell if the energy models are different in H->compute_energy_restricted_emodel
                // requires there to be exactly 2 energy models
				if( is_cross_model(i,j)) {   // If cross model
                   min_en[0] += START_HYBRID_PENALTY;                                            
				}
			}
            //16 Aug 2017 kevin
            //added if is_weakly_closed_flag to make sure there are no base pair going out of region ij, aka is_weakly_closed
            //have it as a flag instead of doing the check here because is_weakly_closed(int i, int j) is a method from the class pseudo_loop which is out of this class
            //so we have to do the check before we call this function which is W_final, a class containing a pseudo_loop object
            if(is_weakly_closed_flag){
                for (auto &energy_model : *energy_models) {
                    energy_model.energy_value = S->compute_energy_restricted_pkonly_emodel (i, j, fres, &energy_model);
                }
                min_en[1] = emodel_energy_function (i, j, energy_models);

                for (auto &energy_model : *energy_models) {
                    energy_model.energy_value = VBI->compute_energy_restricted_pkonly_emodel (i, j, fres, &energy_model);
                }
                min_en[2] = emodel_energy_function (i, j, energy_models);

                for (auto &energy_model : *energy_models) {
                    energy_model.energy_value = VM->compute_energy_restricted_emodel (i, j, fres, &energy_model); // should be left as is, Hosna April 18, 2012
                }
                min_en[3] = emodel_energy_function (i, j, energy_models);

            }
        }
    }

    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = STACK; break;
        case  2: type = INTER; break;
        case  3: type = MULTI; break;
        default: type = NONE;
    }

    if (min_rank > -1)
    {
        if (debug)
            printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
    }

    if (min < INF/2) {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}


// Hosna, April 18, 2012
// pkonly version
void s_energy_matrix::compute_energy_restricted_pkonly (int i, int j, str_features *fres)
// compute the V(i,j) value, if the structure must be restricted
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

	if (fres[i].pair == j && fres[j].pair ==i)
    {
        if (fres[i].pair == i+1)
            min_en[0] = 0;
        else
        {
            if (!exists_restricted (i, j, fres))
                min_en[0] = H->compute_energy_restricted (i, j, fres);

            min_en[1] = S->compute_energy_restricted_pkonly (i, j,fres);
            min_en[2] = VBI->compute_energy_restricted_pkonly (i, j, fres);
            min_en[3] = VM->compute_energy_restricted (i, j, fres); // should be left as is, Hosna April 18, 2012
        }
    }

    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = STACK; break;
        case  2: type = INTER; break;
        case  3: type = MULTI; break;
        default: type = NONE;
    }

    if (min_rank > -1)
    {
        if (debug)
            printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
    }

    if (min < INF/2)
    {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}


void s_energy_matrix::compute_energy_sub (int i, int j)
// suboptimals computation for V(i,j)
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;
    int ij;

    min = INF;
    min_rank = -1;


    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

    if (can_pair (sequence[i], sequence[j]))
    {
        min_en[0] = H->compute_energy (i, j);
        if (i<=j-TURN-1)
        {
            min_en[1] = S->compute_energy (i, j);
            // TODO: uncomment
            if (!ignore_internal)
                min_en[2] = VBI->compute_energy (i, j);
            if (!ignore_multi)
                min_en[3] = VM_sub->compute_energy (i, j);
        }
    }

    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
      {
      case  0: type = HAIRP; break;
      case  1: type = STACK; break;
      case  2: type = INTER; break;
      case  3: type = MULTI; break;
      default: type = NONE;
    }

    if (min < INF/2)
    {
        ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
        //    printf ("V(%d,%d) = %d, type=%c\n", i,j, nodes[ij].energy, nodes[ij].type);
    }
}




void s_energy_matrix::compute_energy_sub_restricted (int i, int j, str_features *fres)
// compute the V(i,j) value - suboptimals and restricted
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

    if (can_pair (sequence[i], sequence[j]))
    {
        min_en[0] = H->compute_energy_restricted (i, j, fres);
        min_en[1] = S->compute_energy (i, j);
        min_en[2] = VBI->compute_energy_restricted (i, j, fres);
        // I don't need restricted for VM_sub because I include dangling ends all the time
        min_en[3] = VM_sub->compute_energy (i, j);
    }

    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = STACK; break;
        case  2: type = INTER; break;
        case  3: type = MULTI; break;
        default: type = NONE;
    }

    if (min_rank > -1)
    {
        if (debug)
            printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
    }

    if (min < INF/2)
    {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}


//kevin 4 oct 2017
void s_energy_matrix::compute_hotspot_energy (int i, int j, int is_stack)
{
    //printf("in compute_hotspot_energy i:%d j:%d\n",i,j);
    PARAMTYPE energy = 0;
    if(is_stack){
        energy = S->compute_energy (i, j);
        //printf("stack: %d\n",energy);
    }else{
        energy = 0;
        // energy = H->compute_energy (i, j);
        //printf("hairpin: %d\n",energy);
    }
        
    //printf ("V(%d,%d) is_stack: %d energy %d\n", i, j, is_stack, energy);
    
    int ij = index[i]+j-i;
    nodes[ij].energy = energy;
    return;
}