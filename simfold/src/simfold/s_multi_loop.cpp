/***************************************************************************
                          s_multi_loop.cpp  -  description
                             -------------------
    begin                : Mon Apr 15 2002
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

// The main class for multi-loop related functions, the mfe case

#include <stdio.h>
//Hosna, March 5, 2012
// malloc.h is not needed on my mac as stdlib.h does the same
//#include <malloc.h>

#include "externs.h"
#include "common.h"
#include "s_multi_loop.h"
#include "simfold.h"
#include "s_energy_matrix.h"
#include "constants.h"

s_multi_loop::s_multi_loop (int *seq, int length)
// The constructor
{
    int i;
    sequence = seq;
    seqlen = length;
    this->V = NULL;

    index = new int[length];    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    int total_length = (length *(length+1))/2;
    index[0] = 0;
    for (i=1; i < length; i++)
        index[i] = index[i-1]+length-i+1;

    WM = new PARAMTYPE [total_length];
    if (WM == NULL) giveup ("Cannot allocate memory", "s_multi_loop");
    for (i=0; i < total_length; i++) WM[i] = INF;
}



s_multi_loop::~s_multi_loop ()
// The destructor
{
    delete [] index;
    delete [] WM;
}


void s_multi_loop::compute_energy_WM (int j)
// compute de MFE of a partial multi-loop closed at (i,j)
{
    int i;
    PARAMTYPE tmp;

    for (i=j-TURN-1; i>=0; i--)
    {
        int ij = index[i]+j-i;
        int iplus1j = index[i+1]+j-i-1;
        int ijminus1 = index[i]+j-1-i;


        tmp = V->get_energy(i,j) +
                   AU_penalty (sequence[i], sequence[j]) +
                   misc.multi_helix_penalty;
        if (tmp < WM[ij])
          {
            WM[ij] = tmp;
          }
        tmp = V->get_energy(i+1,j) +
                  AU_penalty (sequence[i+1], sequence[j]) +
                  dangle_bot [sequence[j]]
                             [sequence[i+1]]
                             [sequence[i]] +
                  misc.multi_helix_penalty +
                  misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i] = -1;
            tmp = tmp - loss (i,i);
        }
        if (tmp < WM[ij])
        {
            WM[ij] = tmp;
        }

        tmp = V->get_energy(i,j-1) +
                  AU_penalty (sequence[i], sequence[j-1]) +
                  dangle_top [sequence [j-1]]
                             [sequence [i]]
                             [sequence [j]] +
                  misc.multi_helix_penalty +
                  misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[j] = -1;
            tmp = tmp - loss (j,j);
        }
        if (tmp < WM[ij])
        {
            WM[ij] = tmp;
        }

        tmp = V->get_energy(i+1,j-1) +
                  AU_penalty (sequence[i+1], sequence[j-1]) +
                  dangle_bot [sequence[j-1]]
                             [sequence[i+1]]
                             [sequence[i]] +
                  dangle_top [sequence [j-1]]
                             [sequence [i+1]]
                             [sequence [j]] +
                  misc.multi_helix_penalty +
                  2*misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i] = -1;
            pred_pairings[j] = -1;
            tmp = tmp - loss (i,i) - loss(j,j);
        }
        if (tmp < WM[ij])
        {
                WM[ij] = tmp;
        }

        tmp = WM[iplus1j] + misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i] = -1;
            tmp = tmp - loss (i,i);
        }
        if (tmp < WM[ij])
        {
            WM[ij] = tmp;
        }

        tmp = WM[ijminus1] + misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[j] = -1;
            tmp = tmp - loss(j,j);
        }
        if (tmp < WM[ij])
        {
            WM[ij] = tmp;
        }

        for (int k=i; k < j; k++)
        {
            int ik = index[i]+k-i;
            int kplus1j = index[k+1]+j-k-1;
            tmp = WM[ik] + WM[kplus1j];
            if (tmp < WM[ij])
              {
                WM[ij] = tmp;
              }
        }
    }
}



PARAMTYPE s_multi_loop::compute_energy (int i, int j)
// compute the MFE of a multi-loop closed at (i,j)
{
    PARAMTYPE min = INF, tmp;
    int k;
    int iplus1k;
    int kplus1jminus1;
    int iplus2k;
    int kplus1jminus2;

    for (k = i+TURN+1; k <= j-TURN-2; k++)
    {
        iplus1k = index[i+1] + k -i-1;
        kplus1jminus1 = index[k+1] + j-1 -k-1;
        iplus2k = index[i+2] + k -i-2;
        kplus1jminus2 = index[k+1] + j-2 -k-1;

        tmp = WM[iplus1k] + WM[kplus1jminus1];
        if (tmp < min)
            min = tmp;

        tmp = WM[iplus2k] + WM[kplus1jminus1] +
              dangle_top [sequence [i]]
              [sequence [j]]
              [sequence [i+1]] +
              misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i+1] = -1;
            tmp = tmp - loss (i+1,i+1);
        }
        if (tmp < min)
            min = tmp;

        tmp = WM[iplus1k] + WM[kplus1jminus2] +
              dangle_bot [sequence[i]]
              [sequence[j]]
              [sequence[j-1]] +
              misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[j-1] = -1;
            tmp = tmp - loss (j-1,j-1);
        }
        if (tmp < min)
            min = tmp;

        tmp = WM[iplus2k] + WM[kplus1jminus2] +
              dangle_top [sequence [i]]
              [sequence [j]]
              [sequence [i+1]] +
              dangle_bot [sequence[i]]
              [sequence[j]]
              [sequence[j-1]] +
              2 * misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i+1] = -1;
            pred_pairings[j-1] = -1;
            tmp = tmp - loss (i+1,i+1) - loss (j-1,j-1);
        }
        if (tmp < min)
          min = tmp;
    }

    min += misc.multi_helix_penalty + misc.multi_offset +
           AU_penalty (sequence[i], sequence[j]);
    return min;
}

void s_multi_loop::compute_energy_WM_restricted (int j, str_features *fres)
// compute de MFE of a partial multi-loop closed at (i,j), the restricted case
{
    int i;
    PARAMTYPE tmp;

    for (i=j-1; i>=0; i--)
    {
        int ij = index[i]+j-i;
        int iplus1j = index[i+1]+j-i-1;
        int ijminus1 = index[i]+j-1-i;

        tmp = V->get_energy(i,j) +
              AU_penalty (sequence[i], sequence[j]) +
              misc.multi_helix_penalty;

        if (tmp < WM[ij]) {
            WM[ij] = tmp;
        }

        if (fres[i].pair <= -1)
        {
            tmp = V->get_energy(i+1,j) +
                    AU_penalty (sequence[i+1], sequence[j]) +
                    dangle_bot [sequence[j]]
                                [sequence[i+1]]
                                [sequence[i]] +
                    misc.multi_helix_penalty +
                    misc.multi_free_base_penalty;
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }
        if (fres[j].pair <= -1)
        {
            tmp = V->get_energy(i,j-1) +
                    AU_penalty (sequence[i], sequence[j-1]) +
                    dangle_top [sequence [j-1]]
                                [sequence [i]]
                                [sequence [j]] +
                    misc.multi_helix_penalty +
                    misc.multi_free_base_penalty;
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }

        if (fres[i].pair <= -1 && fres[j].pair <= -1)
        {
            tmp = V->get_energy(i+1,j-1) +
                    AU_penalty (sequence[i+1], sequence[j-1]) +
                    dangle_bot [sequence[j-1]]
                                [sequence[i+1]]
                                [sequence[i]] +
                    dangle_top [sequence [j-1]]
                                [sequence [i+1]]
                                [sequence [j]] +
                    misc.multi_helix_penalty +
                    2*misc.multi_free_base_penalty;
            if (tmp < WM[ij])
            {
                    WM[ij] = tmp;
            }
        }

        if (fres[i].pair <= -1)
        {
            tmp = WM[iplus1j] + misc.multi_free_base_penalty;
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }

        if (fres[j].pair <= -1)
        {
            tmp = WM[ijminus1] + misc.multi_free_base_penalty;
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }

        for (int k=i; k < j; k++)
        {
            int ik = index[i]+k-i;
            int kplus1j = index[k+1]+j-k-1;
            tmp = WM[ik] + WM[kplus1j];
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }
    }
}

void s_multi_loop::compute_energy_WM_restricted_pmo (int j, str_features *fres)
// compute de MFE of a partial multi-loop closed at (i,j), the restricted case
{
    int i;
    PARAMTYPE tmp;
	double pmo_percentage, rna_percentage;

    for (i=j-1; i>=0; i--)
    {
		get_pmo_usage_percentages(i, j, &pmo_percentage, &rna_percentage);
		if (pmo_percentage == 0 && rna_percentage == 0)
			continue;

        int ij = index[i]+j-i;
        int iplus1j = index[i+1]+j-i-1;
        int ijminus1 = index[i]+j-1-i;

        tmp = V->get_energy(i,j) +
              (PARAMTYPE) (pmo_percentage*AU_penalty_pmo(sequence[i], sequence[j]) + rna_percentage*AU_penalty(sequence[i], sequence[j])) +
              (PARAMTYPE) (pmo_percentage*misc_pmo.multi_helix_penalty + rna_percentage*misc.multi_helix_penalty);

        if (tmp < WM[ij]) {
            WM[ij] = tmp;
        }

        if (fres[i].pair <= -1)
        {
            tmp = V->get_energy(i+1,j) +
                  (PARAMTYPE) (pmo_percentage*AU_penalty_pmo(sequence[i+1], sequence[j]) + rna_percentage*AU_penalty(sequence[i+1], sequence[j])) +
                  (PARAMTYPE) (pmo_percentage*dangle_bot_pmo[sequence[j]][sequence[i+1]][sequence[i]] + rna_percentage*dangle_bot[sequence[j]][sequence[i+1]][sequence[i]]) +
                  (PARAMTYPE) (pmo_percentage*misc_pmo.multi_helix_penalty + rna_percentage*misc.multi_helix_penalty); +
                  (PARAMTYPE) (pmo_percentage*misc_pmo.multi_free_base_penalty + rna_percentage*misc.multi_free_base_penalty);

			if (tmp < WM[ij]) {
                WM[ij] = tmp;
            }
        }
        if (fres[j].pair <= -1)
        {
            tmp = V->get_energy(i,j-1) +
                  (PARAMTYPE) (pmo_percentage*AU_penalty_pmo(sequence[i], sequence[j-1]) + rna_percentage*AU_penalty(sequence[i], sequence[j-1])) +
                  (PARAMTYPE) (pmo_percentage*dangle_top_pmo[sequence [j-1]][sequence [i]][sequence [j]] + rna_percentage*dangle_top[sequence [j-1]][sequence [i]][sequence [j]]) +
                  (PARAMTYPE) (pmo_percentage*misc_pmo.multi_helix_penalty + rna_percentage*misc.multi_helix_penalty); +
                  (PARAMTYPE) (pmo_percentage*misc_pmo.multi_free_base_penalty + rna_percentage*misc.multi_free_base_penalty);

            if (tmp < WM[ij]) {
                WM[ij] = tmp;
            }
        }

        if (fres[i].pair <= -1 && fres[j].pair <= -1)
        {
            tmp = V->get_energy(i+1,j-1) +
                  (PARAMTYPE) (pmo_percentage*AU_penalty_pmo(sequence[i+1], sequence[j-1]) + rna_percentage*AU_penalty(sequence[i+1], sequence[j-1])) +
                  (PARAMTYPE) (pmo_percentage*dangle_bot_pmo[sequence[j-1]][sequence[i+1]][sequence[i]] + rna_percentage*dangle_bot[sequence[j-1]][sequence[i+1]][sequence[i]]) +
                  (PARAMTYPE) (pmo_percentage*dangle_top_pmo[sequence [j-1]][sequence [i+1]][sequence [j]] + rna_percentage*dangle_top[sequence [j-1]][sequence [i+1]][sequence [j]]) +
                  (PARAMTYPE) (pmo_percentage*misc_pmo.multi_helix_penalty + rna_percentage*misc.multi_helix_penalty); +
                  (PARAMTYPE) (pmo_percentage*2*misc_pmo.multi_free_base_penalty + rna_percentage*2*misc.multi_free_base_penalty);

            if (tmp < WM[ij]) {
            	WM[ij] = tmp;
            }
        }

        if (fres[i].pair <= -1)
        {
            tmp = WM[iplus1j] + (PARAMTYPE) (pmo_percentage*misc_pmo.multi_free_base_penalty + rna_percentage*misc.multi_free_base_penalty);

            if (tmp < WM[ij]) {
                WM[ij] = tmp;
            }
        }

        if (fres[j].pair <= -1)
        {
            tmp = WM[ijminus1] + (PARAMTYPE) (pmo_percentage*misc_pmo.multi_free_base_penalty + rna_percentage*misc.multi_free_base_penalty);

            if (tmp < WM[ij]) {
                WM[ij] = tmp;
            }
        }

        for (int k=i; k < j; k++)
        {
            int ik = index[i]+k-i;
            int kplus1j = index[k+1]+j-k-1;
            tmp = WM[ik] + WM[kplus1j];
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }
    }
}

PARAMTYPE s_multi_loop::compute_energy_restricted (int i, int j, str_features *fres)
// compute the MFE of a multi-loop closed at (i,j), the restricted case
{
    PARAMTYPE min = INF, tmp;
    int k;
    int iplus1k;
    int kplus1jminus1;
    int iplus2k;
    int kplus1jminus2;

	if (sequence[i] == 4 || sequence[j] == 4)
		return 0;

    // May 16, 2007: Replaced this for loop, because we may have very short restricted branches
    //for (k = i+TURN+1; k <= j-TURN-2; k++)
    for (k = i+2; k <= j-3; k++)
    {
        iplus1k = index[i+1] + k -i-1;
        kplus1jminus1 = index[k+1] + j-1 -k-1;
        iplus2k = index[i+2] + k -i-2;
        kplus1jminus2 = index[k+1] + j-2 -k-1;

        tmp = WM[iplus1k] + WM[kplus1jminus1];
        if (tmp < min)
            min = tmp;

        if (fres[i+1].pair <= -1)
        {
            tmp = WM[iplus2k] + WM[kplus1jminus1] +
                dangle_top [sequence [i]]
                [sequence [j]]
                [sequence [i+1]] +
                misc.multi_free_base_penalty;
            if (tmp < min)
                min = tmp;
        }
        if (fres[j-1].pair <= -1)
        {
            tmp = WM[iplus1k] + WM[kplus1jminus2] +
                dangle_bot [sequence[i]]
                [sequence[j]]
                [sequence[j-1]] +
                misc.multi_free_base_penalty;
            if (tmp < min)
                min = tmp;
        }
        if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
        {
            tmp = WM[iplus2k] + WM[kplus1jminus2] +
                dangle_top [sequence [i]]
                [sequence [j]]
                [sequence [i+1]] +
                dangle_bot [sequence[i]]
                [sequence[j]]
                [sequence[j-1]] +
                2 * misc.multi_free_base_penalty;
            if (tmp < min)
            min = tmp;
        }
    }

    min += misc.multi_helix_penalty + misc.multi_offset +
           AU_penalty (sequence[i], sequence[j]);
    return min;
}

//AP
PARAMTYPE s_multi_loop::compute_energy_restricted_emodel (int i, int j, str_features *fres, energy_model *model)
// compute the MFE of a multi-loop closed at (i,j), the restricted case
{
    PARAMTYPE min = INF, tmp;
    int k;
    int iplus1k;
    int kplus1jminus1;
    int iplus2k;
    int kplus1jminus2;

/*
	//AP
	if (sequence[i] == X || sequence[j] == X)
		return 0; // Why 0 not INF?
	if (sequence[i+1] == X || sequence[j+1] == X || sequence[i-1] == X || sequence[j-1] == X)
		return INF;
*/

    // Ian Wark and Kevin July 20 2017
    // if i or j is linker (X), cannot be anything
    // Uses i+1 or j-1, so those also cannot be X
    if (sequence[i] == X || sequence[j] == X)
        return INF;
    
    //handles i+1 and j-1 are both X, so we have to force haipin to be picked
    //also prevents i and j to cross each other
    if(sequence[i+1] == X && sequence[j-1] == X){
        return INF;
    }

    //if i+1 is a X, then we have to shift the beginning of k to the index after the X
    //shift i 
    if(sequence[i+1] == X){
        i++;
        while(sequence[i] == X){
            i++;
        }    
    } 

    //if j-1 is a X, then we have to shift the end of k to the index before the X
    //shift j 
    if(sequence[j-1] == X){
        j--;
        while(sequence[j] == X){
            j--;
        }    
    }

	// May 16, 2007: Replaced this for loop, because we may have very short restricted branches
    //for (k = i+TURN+1; k <= j-TURN-2; k++)
    for (k = i+2; k <= j-3; k++)
    {
        iplus1k = index[i+1] + k -i-1;
        kplus1jminus1 = index[k+1] + j-1 -k-1;
        iplus2k = index[i+2] + k -i-2;
        kplus1jminus2 = index[k+1] + j-2 -k-1;

        tmp = WM[iplus1k] + WM[kplus1jminus1];
        if (tmp < min)
            min = tmp;

        if (fres[i+1].pair <= -1)
        {
            tmp = WM[iplus2k] + WM[kplus1jminus1] +
                model->dangle_top [sequence [i]]
								[sequence [j]]
								[sequence [i+1]] +
                model->misc.multi_free_base_penalty; 
            if (tmp < min)
                min = tmp;
        }
        if (fres[j-1].pair <= -1)
        {
            tmp = WM[iplus1k] + WM[kplus1jminus2] +
                model->dangle_bot [sequence[i]]
								[sequence[j]]
								[sequence[j-1]] +
                model->misc.multi_free_base_penalty;
            if (tmp < min)
                min = tmp;
        }
        if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
        {
            tmp = WM[iplus2k] + WM[kplus1jminus2] +
                model->dangle_top [sequence [i]]
								[sequence [j]]
								[sequence [i+1]] +
                model->dangle_bot [sequence[i]]
								[sequence[j]]
								[sequence[j-1]] +
                2 * model->misc.multi_free_base_penalty;
            if (tmp < min)
            	min = tmp;
        }
    }

    min += model->misc.multi_helix_penalty + model->misc.multi_offset +
           AU_penalty_emodel (sequence[i], sequence[j], model);
    return min;
}

//AP
void s_multi_loop::compute_energy_WM_restricted_emodel (int j, str_features *fres, std::vector<energy_model> *energy_models) {
// compute de MFE of a partial multi-loop closed at (i,j), the restricted case
    int i;
    PARAMTYPE tmp;
	energy_model *model;
    //Aug 14 kevin and Ian and Mahyar and Hosna
    //used for handling X cases
    int temp_ij = -1;
    int temp_i = 0;
    int temp_j = 0;

    for (i=j-1; i>=0; i--) {
        int ij = index[i]+j-i;
        int iplus1j = index[i+1]+j-i-1;
        int ijminus1 = index[i]+j-1-i;

        //Aug 14 2017 kevin and Ian and Mayhar and Hosna
        //WM[ij] = INF if i,j both X or ij cross each other
        if(sequence[i] == X && sequence[j] == X){
            WM[ij] = INF;
            continue;
        }

        //Aug 14 2017 kevin and Ian and Mayhar and Hosna
        //if i is X, shift i till i is not a X and use that value
        //if temp_i,j cross each other, WM[ij] = INF
        if (sequence[i] == X){
            temp_i = i;
            while(sequence[temp_i] == X){
                temp_i++;
            }
            if(temp_i < j){
                temp_ij = index[temp_i]+j-temp_i;
                WM[ij] = WM[temp_ij]; 
            }else{
                WM[ij] = INF;
            }
            continue;
        }

        //Aug 14 2017 kevin and Ian and Mayhar and Hosna
        //if j is X, shift j till j is not a X and use that value
        //if i,temp_j cross each other, WM[ij] = INF
        if (sequence[j] == X) {
            temp_j = j;
            while(sequence[temp_j] == X){
                temp_j--;
            }
            if(i < temp_j){
                temp_ij = index[i]+temp_j-i;
                WM[ij] = WM[temp_ij];
            }else{
                WM[ij] = INF;
            }
            continue;
	    }

        //normal calculations without X

		//AP
		for (auto &energy_model : *energy_models) {
			model = &energy_model;
			model->energy_value = WM[ij];

		    tmp = V->get_energy(i,j) +
		          AU_penalty_emodel (sequence[i], sequence[j], model) +
		          model->misc.multi_helix_penalty;
		    if (tmp < model->energy_value) {
		        model->energy_value = tmp;
		    }

            if (fres[i].pair <= -1 && fres[j].pair <= -1)
		    {
		        tmp = V->get_energy(i+1,j-1) +
		                AU_penalty_emodel (sequence[i+1], sequence[j-1], model) +
		                model->misc.multi_helix_penalty +
		                2*model->misc.multi_free_base_penalty;
                //Aug 14 2017 kevin and Mahyar
                //modiefied the formula such that we only add dangle_bot,dangle_top when i,j,i+1,j-1 are not X to avoid seg fault
                if(sequence[i] != X && sequence[j] != X && sequence[i+1] != X && sequence[j-1] != X){
                    tmp += model->dangle_bot [sequence[j-1]]
		                            	[sequence[i+1]]
		                            	[sequence[i]] +
		                model->dangle_top [sequence [j-1]]
		                            	[sequence [i+1]]
		                            	[sequence [j]];
                }
		        if (tmp < model->energy_value)
		        {
		                model->energy_value = tmp;
		        }
		    }

		    if (fres[i].pair <= -1)
		    {
                tmp = V->get_energy(i+1,j) +
                        AU_penalty_emodel (sequence[i+1], sequence[j], model) +
                        model->misc.multi_helix_penalty +
                        model->misc.multi_free_base_penalty;
                //Aug 14 2017 kevin and Mahyar
                //modiefied the formula such that we only add dangle_bot when i,j,i+1 are not X to avoid seg fault
                if(sequence[i] != X && sequence[j] != X && sequence[i+1] != X){
                    tmp += model->dangle_bot [sequence[j]]
                                        [sequence[i+1]]
                                        [sequence[i]];
                }
                if (tmp < model->energy_value)
                {
                    model->energy_value = tmp;
                }
                
                
		    }
		    if (fres[j].pair <= -1)
		    {
		        tmp = V->get_energy(i,j-1) +
		                AU_penalty_emodel (sequence[i], sequence[j-1], model) +
		                model->misc.multi_helix_penalty +
		                model->misc.multi_free_base_penalty;
                //Aug 14 2017 kevin and Mahyar
                //modiefied the formula such that we only add dangle_top when i,j,j-1 are not X to avoid seg fault
                if(sequence[i] != X && sequence[j] != X && sequence[j-1] != X){
                    tmp += model->dangle_top [sequence [j-1]]
		                            	[sequence [i]]
		                            	[sequence [j]];
                }
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		        }
		    }

		    if (fres[i].pair <= -1)
		    {
		        tmp = WM[iplus1j] + model->misc.multi_free_base_penalty;
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		        }
		    }

		    if (fres[j].pair <= -1)
		    {
		        tmp = WM[ijminus1] + model->misc.multi_free_base_penalty;
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		        }
		    }

		    for (int k=i; k < j; k++)
		    {
                if (sequence[k] == X)
                    continue;

		        int ik = index[i]+k-i;
		        int kplus1j = index[k+1]+j-k-1;
		        tmp = WM[ik] + WM[kplus1j];
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		        }
		    }
		}

		WM[ij] = emodel_energy_function (i, j, energy_models);
    }
}

//AP
void s_multi_loop::compute_energy_WM_restricted_pkonly_emodel (int j, str_features *fres, std::vector<energy_model> *energy_models)
// compute the MFE of a partial multi-loop closed at (i,j), the restricted case when i and j are already paired
{
    int i;
    PARAMTYPE tmp = INF;
	energy_model *model;
    //Aug 14 kevin and Ian and Mahyar and Hosna
    //used for handling X cases
    int temp_ij = -1;
    int temp_i = 0;
    int temp_j = 0;

    for (i=j-1; i>=0; i--) {
        int ij = index[i]+j-i;
        int iplus1j = index[i+1]+j-i-1;
        int ijminus1 = index[i]+j-1-i;
        
        //Aug 14 2017 kevin and Ian and Mayhar and Hosna
        //WM[ij] = INF if i,j both X or ij cross each other
        if(sequence[i] == X && sequence[j] == X){
            WM[ij] = INF;
            continue;
        }

        //Aug 14 2017 kevin and Ian and Mayhar and Hosna
        //if i is X, shift i till i is not a X and use that value
        //if temp_i,j cross each other, WM[ij] = INF
        if (sequence[i] == X){
            temp_i = i;
            while(sequence[temp_i] == X){
                temp_i++;
            }
            if(temp_i < j){
                temp_ij = index[temp_i]+j-temp_i;
                WM[ij] = WM[temp_ij]; 
            }else{
                WM[ij] = INF;
            }
            continue;
        }

        //Aug 14 2017 kevin and Ian and Mayhar and Hosna
        //if j is X, shift j till j is not a X and use that value
        //if i,temp_j cross each other, WM[ij] = INF
        if (sequence[j] == X) {
            temp_j = j;
            while(sequence[temp_j] == X){
                temp_j--;
            }
            if(i < temp_j){
                temp_ij = index[i]+temp_j-i;
                WM[ij] = WM[temp_ij];
            }else{
                WM[ij] = INF;
            }
            continue;
	    }

        //normal calculations without X

		//AP
		for (auto &energy_model : *energy_models) {
			model = &energy_model;
			model->energy_value = WM[ij];

			if (fres[i].pair == j && fres[j].pair == i){
				tmp = V->get_energy(i,j) +
				AU_penalty_emodel (sequence[i], sequence[j], model) +
				model->misc.multi_helix_penalty;
			}

		    if (tmp < model->energy_value)
			{
		        model->energy_value = tmp;
			}

		    if (fres[i].pair <= -1 && fres[i+1].pair == j && fres[j].pair == i+1)
		    {
		        tmp = V->get_energy(i+1,j) +
				AU_penalty_emodel (sequence[i+1], sequence[j], model) +
				model->misc.multi_helix_penalty +
				model->misc.multi_free_base_penalty;
                //Aug 14 2017 kevin and Mahyar
                //modiefied the formula such that we only add dangle_bot when i,j,i+1 are not X to avoid seg fault
                if(sequence[i] != X && sequence[i+1] != X && sequence[j] != X){
                    tmp += model->dangle_bot [sequence[j]]
								[sequence[i+1]]
								[sequence[i]];
                }
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		        }
		    }
		    if (fres[j].pair <= -1 && fres[i].pair == j-1 && fres[j-1].pair == i)
		    {
		        tmp = V->get_energy(i,j-1) +
				AU_penalty_emodel (sequence[i], sequence[j-1], model) +
				model->misc.multi_helix_penalty +
				model->misc.multi_free_base_penalty;
                //Aug 14 2017 kevin and Mahyar
                //modiefied the formula such that we only add dangle_top when i,j,j-1 are not X to avoid seg fault
                if(sequence[i] != X && sequence[j] != X && sequence[j-1] != X){
                    model->dangle_top [sequence [j-1]]
								[sequence [i]]
								[sequence [j]];
                }
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		        }
		    }

		    if (fres[i].pair <= -1 && fres[j].pair <= -1 && fres[i+1].pair == j-1 && fres[j-1].pair == i+1)
		    {
		        tmp = V->get_energy(i+1,j-1) +
				AU_penalty_emodel (sequence[i+1], sequence[j-1], model) +
				model->misc.multi_helix_penalty +
				2*model->misc.multi_free_base_penalty;
                //Aug 14 2017 kevin and Mahyar
                //modiefied the formula such that we only add dangle_bot,dangle_top when i,j,i+1,j-1 are not X to avoid seg fault
                if(sequence[i] != X && sequence[j] != X && sequence[i+1] != X && sequence[j-1] != X){
                    tmp += model->dangle_bot [sequence[j-1]]
								[sequence[i+1]]
								[sequence[i]] +
				            model->dangle_top [sequence [j-1]]
								[sequence [i+1]]
								[sequence [j]];
                }
		        if (tmp < model->energy_value)
		        {
					model->energy_value = tmp;
		        }
		    }

		    if (fres[i].pair <= -1)
		    {
		        tmp = WM[iplus1j] + model->misc.multi_free_base_penalty;
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		        }
		    }

		    if (fres[j].pair <= -1)
		    {
		        tmp = WM[ijminus1] + model->misc.multi_free_base_penalty;
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		        }
		    }

		    for (int k=i; k < j; k++)
		    {
                if(sequence[k] == X){
                    continue;
                }
		        int ik = index[i]+k-i;
		        int kplus1j = index[k+1]+j-k-1;
		        tmp = WM[ik] + WM[kplus1j];
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		        }
		    }
		}

		WM[ij] = emodel_energy_function (i, j, energy_models);
    }
}

PARAMTYPE s_multi_loop::compute_energy_restricted_pmo (int i, int j, str_features *fres)
// compute the MFE of a multi-loop closed at (i,j), the restricted case
{
    PARAMTYPE min = INF, tmp;
    int k;
    int iplus1k;
    int kplus1jminus1;
    int iplus2k;
    int kplus1jminus2;
	double pmo_percentage, rna_percentage;

	get_pmo_usage_percentages(i, j, &pmo_percentage, &rna_percentage);
	if (pmo_percentage == 0 && rna_percentage == 0)
		return 0;

    // May 16, 2007: Replaced this for loop, because we may have very short restricted branches
    //for (k = i+TURN+1; k <= j-TURN-2; k++)
    for (k = i+2; k <= j-3; k++)
    {
        iplus1k = index[i+1] + k -i-1;
        kplus1jminus1 = index[k+1] + j-1 -k-1;
        iplus2k = index[i+2] + k -i-2;
        kplus1jminus2 = index[k+1] + j-2 -k-1;

        tmp = WM[iplus1k] + WM[kplus1jminus1];
        if (tmp < min)
            min = tmp;

        if (fres[i+1].pair <= -1)
        {
            tmp = WM[iplus2k] +
				WM[kplus1jminus1] +
                (PARAMTYPE) (pmo_percentage*dangle_top_pmo[sequence[i]]
														[sequence[j]]
														[sequence[i+1]] +
								rna_percentage*dangle_top[sequence[i]]
														[sequence[j]]
														[sequence[i+1]]) +
                (PARAMTYPE) (pmo_percentage*misc_pmo.multi_free_base_penalty + rna_percentage*misc.multi_free_base_penalty);

			if (tmp < min)
                min = tmp;
        }
        if (fres[j-1].pair <= -1)
        {
            tmp = WM[iplus1k] +
				WM[kplus1jminus2] +
                (PARAMTYPE) (pmo_percentage*dangle_bot_pmo[sequence[i]]
														[sequence[j]]
														[sequence[j-1]] +
								rna_percentage*dangle_bot[sequence[i]]
														[sequence[j]]
														[sequence[j-1]]) +
                (PARAMTYPE) (pmo_percentage*misc_pmo.multi_free_base_penalty + rna_percentage*misc.multi_free_base_penalty);

			if (tmp < min)
                min = tmp;
        }
        if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
        {
            tmp = WM[iplus2k] +
				WM[kplus1jminus2] +
                (PARAMTYPE) (pmo_percentage*dangle_top_pmo[sequence[i]]
														[sequence[j]]
														[sequence[i+1]] +
								rna_percentage*dangle_top[sequence[i]]
														[sequence[j]]
														[sequence[i+1]]) +
                (PARAMTYPE) (pmo_percentage*dangle_bot_pmo[sequence[i]]
														[sequence[j]]
														[sequence[j-1]] +
								rna_percentage*dangle_bot[sequence[i]]
														[sequence[j]]
														[sequence[j-1]]) +
                (PARAMTYPE) (pmo_percentage*2*misc_pmo.multi_free_base_penalty + rna_percentage*2*misc.multi_free_base_penalty);

			if (tmp < min)
            	min = tmp;
        }
    }

    min += (PARAMTYPE) (pmo_percentage*(misc_pmo.multi_helix_penalty +
										misc_pmo.multi_offset +
										AU_penalty_pmo(sequence[i], sequence[j])) +
						rna_percentage*(misc.multi_helix_penalty +
										misc.multi_offset +
										AU_penalty(sequence[i], sequence[j])));

	return min;
}

// added April 18, 2012
void s_multi_loop::compute_energy_WM_restricted_pkonly (int j, str_features *fres)
// compute the MFE of a partial multi-loop closed at (i,j), the restricted case when i and j are already paired
{
    int i;
    PARAMTYPE tmp = INF;

    for (i=j-1; i>=0; i--)
    {
        int ij = index[i]+j-i;
        int iplus1j = index[i+1]+j-i-1;
        int ijminus1 = index[i]+j-1-i;

		if (fres[i].pair == j && fres[j].pair == i){
			tmp = V->get_energy(i,j) +
			AU_penalty (sequence[i], sequence[j]) +
			misc.multi_helix_penalty;
		}

        if (tmp < WM[ij])
		{
            WM[ij] = tmp;
		}

        if (fres[i].pair <= -1 && fres[i+1].pair == j && fres[j].pair == i+1)
        {
            tmp = V->get_energy(i+1,j) +
			AU_penalty (sequence[i+1], sequence[j]) +
			dangle_bot [sequence[j]]
			[sequence[i+1]]
			[sequence[i]] +
			misc.multi_helix_penalty +
			misc.multi_free_base_penalty;
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }
        if (fres[j].pair <= -1 && fres[i].pair == j-1 && fres[j-1].pair == i)
        {
            tmp = V->get_energy(i,j-1) +
			AU_penalty (sequence[i], sequence[j-1]) +
			dangle_top [sequence [j-1]]
			[sequence [i]]
			[sequence [j]] +
			misc.multi_helix_penalty +
			misc.multi_free_base_penalty;
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }

        if (fres[i].pair <= -1 && fres[j].pair <= -1 && fres[i+1].pair == j-1 && fres[j-1].pair == i+1)
        {
            tmp = V->get_energy(i+1,j-1) +
			AU_penalty (sequence[i+1], sequence[j-1]) +
			dangle_bot [sequence[j-1]]
			[sequence[i+1]]
			[sequence[i]] +
			dangle_top [sequence [j-1]]
			[sequence [i+1]]
			[sequence [j]] +
			misc.multi_helix_penalty +
			2*misc.multi_free_base_penalty;
            if (tmp < WM[ij])
            {
				WM[ij] = tmp;
            }
        }

        if (fres[i].pair <= -1)
        {
            tmp = WM[iplus1j] + misc.multi_free_base_penalty;
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }

        if (fres[j].pair <= -1)
        {
            tmp = WM[ijminus1] + misc.multi_free_base_penalty;
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }

        for (int k=i; k < j; k++)
        {
            int ik = index[i]+k-i;
            int kplus1j = index[k+1]+j-k-1;
            tmp = WM[ik] + WM[kplus1j];
            if (tmp < WM[ij])
            {
                WM[ij] = tmp;
            }
        }
    }
}

//TODO: PK
void s_multi_loop::compute_energy_WM_restricted_pkonly_pmo (int j, str_features *fres)
// compute the MFE of a partial multi-loop closed at (i,j), the restricted case when i and j are already paired
{
    int i;
	double pmo_percentage, rna_percentage;
    PARAMTYPE tmp = INF;

    for (i=j-1; i>=0; i--) {
        int ij = index[i]+j-i;
        int iplus1j = index[i+1]+j-i-1;
        int ijminus1 = index[i]+j-1-i;

		get_pmo_usage_percentages(i, j, &pmo_percentage, &rna_percentage);
		if (pmo_percentage == 0 && rna_percentage == 0)
			continue;

		if (fres[i].pair == j && fres[j].pair == i) {
			tmp = V->get_energy(i,j) +
			(PARAMTYPE) (pmo_percentage*AU_penalty_pmo(sequence[i], sequence[j]) + rna_percentage*AU_penalty(sequence[i], sequence[j])) +
			(PARAMTYPE) (pmo_percentage*misc_pmo.multi_helix_penalty + rna_percentage*misc.multi_helix_penalty);
		}

        if (tmp < WM[ij]) {
            WM[ij] = tmp;
		}

        if (fres[i].pair <= -1 && fres[i+1].pair == j && fres[j].pair == i+1) {
            tmp = V->get_energy(i+1,j) +
			(PARAMTYPE) (pmo_percentage*AU_penalty_pmo(sequence[i+1], sequence[j]) + rna_percentage*AU_penalty(sequence[i+1], sequence[j])) +
			(PARAMTYPE) (pmo_percentage*dangle_bot_pmo[sequence[j]]
													[sequence[i+1]]
													[sequence[i]] +
							rna_percentage*dangle_bot[sequence[j]]
													[sequence[i+1]]
													[sequence[i]]) +
			(PARAMTYPE) (pmo_percentage*misc_pmo.multi_helix_penalty + rna_percentage*misc.multi_helix_penalty) +
			(PARAMTYPE) (pmo_percentage*misc_pmo.multi_free_base_penalty + rna_percentage*misc.multi_free_base_penalty);

            if (tmp < WM[ij]) {
                WM[ij] = tmp;
            }
        }

        if (fres[j].pair <= -1 && fres[i].pair == j-1 && fres[j-1].pair == i) {
            tmp = V->get_energy(i,j-1) +
			(PARAMTYPE) (pmo_percentage*AU_penalty_pmo(sequence[i], sequence[j-1]) + rna_percentage*AU_penalty(sequence[i], sequence[j-1])) +
			(PARAMTYPE) (pmo_percentage*dangle_top_pmo[sequence [j-1]]
													[sequence [i]]
													[sequence [j]] +
							rna_percentage*dangle_top[sequence [j-1]]
													[sequence [i]]
													[sequence [j]]) +
			(PARAMTYPE) (pmo_percentage*misc_pmo.multi_helix_penalty + rna_percentage*misc.multi_helix_penalty) +
			(PARAMTYPE) (pmo_percentage*misc_pmo.multi_free_base_penalty + rna_percentage*misc.multi_free_base_penalty);

			if (tmp < WM[ij]) {
                WM[ij] = tmp;
            }
        }

        if (fres[i].pair <= -1 && fres[j].pair <= -1 && fres[i+1].pair == j-1 && fres[j-1].pair == i+1) {
            tmp = V->get_energy(i+1,j-1) +
			(PARAMTYPE) (pmo_percentage*AU_penalty_pmo(sequence[i+1], sequence[j-1]) + rna_percentage*AU_penalty(sequence[i+1], sequence[j-1])) +
			(PARAMTYPE) (pmo_percentage*dangle_bot_pmo[sequence[j-1]]
													[sequence[i+1]]
													[sequence[i]] +
							rna_percentage*dangle_bot[sequence[j-1]]
													[sequence[i+1]]
													[sequence[i]]) +

			(PARAMTYPE) (pmo_percentage*dangle_top_pmo[sequence [j-1]]
													[sequence [i+1]]
													[sequence [j]] +
							rna_percentage*dangle_top[sequence [j-1]]
													[sequence [i+1]]
													[sequence [j]]) +
			(PARAMTYPE) (pmo_percentage*misc_pmo.multi_helix_penalty + rna_percentage*misc.multi_helix_penalty) +
			(PARAMTYPE) (pmo_percentage*2*misc_pmo.multi_free_base_penalty + rna_percentage*2*misc.multi_free_base_penalty);

			if (tmp < WM[ij])  {
				WM[ij] = tmp;
            }
        }

        if (fres[i].pair <= -1) {
            tmp = WM[iplus1j] + (PARAMTYPE) (pmo_percentage*misc_pmo.multi_free_base_penalty + rna_percentage*misc.multi_free_base_penalty);

            if (tmp < WM[ij]) {
                WM[ij] = tmp;
            }
        }

        if (fres[j].pair <= -1) {
            tmp = WM[ijminus1] + (PARAMTYPE) (pmo_percentage*misc_pmo.multi_free_base_penalty + rna_percentage*misc.multi_free_base_penalty);

            if (tmp < WM[ij]) {
                WM[ij] = tmp;
            }
        }

        for (int k=i; k < j; k++) {
            int ik = index[i]+k-i;
            int kplus1j = index[k+1]+j-k-1;
            tmp = WM[ik] + WM[kplus1j];

            if (tmp < WM[ij]) {
                WM[ij] = tmp;
            }
        }
    }
}
