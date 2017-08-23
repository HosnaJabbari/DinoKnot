
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "pseudo_loop.h"
#include "V_final.h"
#include "W_final.h"
#include "constants.h"
#include "h_struct.h"
#include "externs.h"
#include "h_externs.h"
#include "h_common.h"


//kevin
#include "simfold.h"
#include "params.h"
#include "structs.h"
#include "common.h"

#define debug 0



// Hosna June 20th, 2007
// calls the constructor for s_min_folding
// to create all the matrixes required for simfold
// and then calls allocate_space in here to allocate
// space for WMB and V_final
W_final::W_final(char *seq, char *res):s_min_folding(seq,res)
{
	this->nb_nucleotides = strlen(seq);

	// Ian Wark July 21 2017
	// we don't need this, this is done in s_min_folding
	/*
	this->int_sequence = new int[this->nb_nucleotides];
	if (int_sequence == NULL) giveup ("Cannot allocate memory", "W_final");
	int i;
    for (i=0; i < this->nb_nucleotides; i++) int_sequence[i] = nuc_to_int(seq[i]);
	*/

	space_allocation();
}

W_final::W_final(char *seq, char *res, std::vector<energy_model> *energy_models):s_min_folding(seq,res,energy_models)
{
	this->nb_nucleotides = strlen(seq);

	// Ian Wark July 21 2017
	// we don't need this, this is done in s_min_folding
	/*
	this->int_sequence = new int[this->nb_nucleotides];
	if (int_sequence == NULL) giveup ("Cannot allocate memory", "W_final");
	int i;
    for (i=0; i < this->nb_nucleotides; i++) int_sequence[i] = nuc_to_int(seq[i]);
    */

	space_allocation();
}

W_final::~W_final()
{
	delete vm;
	delete v;
	delete WMB;

	// Ian Wark July 21 2017
	// we don't need this, this is done in s_min_folding
	//delete [] int_sequence;
}

// Hosna June 20th, 2007
// allocates space for WMB object and V_final
void W_final::space_allocation(){

	// Hosna June 20th, 2007
	vm = new VM_final(this->int_sequence,this->nb_nucleotides);
	if (vm == NULL) giveup ("Cannot allocate memory", "W_final");
	if (debug){
		printf("nb_nucleotides = %d \n",this->nb_nucleotides);
	}

	// Hosna June 20th, 2007
	// I don't think we need the following line
	//vm->set_energy_matrix(s_min_folding::V);

	// Hosna June 20th, 2007
	v = new V_final(nb_nucleotides);
	if (v == NULL) giveup ("Cannot allocate memory", "W_final");
	//s_min_folding::V, s_min_folding::H, s_min_folding::S, s_min_folding::VBI, vm);
	v->setloops(this->V,vm);

	// Hosna: June 20th 2007
    WMB = new pseudo_loop (sequence,restricted,v,this->H,this->S,this->VBI,vm);
    if (WMB == NULL) giveup ("Cannot allocate memory", "W_final");

    // Hosna: June 20th 2007
    vm->set_V_matrix(v);
    vm->set_WMB_matrix(WMB);


}


double W_final::hfold_pkonly(){

	double energy;
    int i, j;

	h_str_features *h_fres;
    if ((h_fres = new h_str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "h_str_features");
    // detect the structure features
    detect_h_structure_features (restricted, h_fres);
    WMB->set_features(h_fres);
    WMB->initialize();

    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    detect_structure_features (restricted, fres);



    // Hosna: June 28, 2007
    // set the features for checking
    v->set_features(fres);

    // Hosna: July 2nd, 2007
    // set the VM matrix for VM_final
    vm->set_VM_matrix(VM);


	// TODO:
	// I think I shoud fill simfold tables here, before filling the HFold tables
	// Hosna, March 8, 2012

	// 1) fill all th ematrices simfold needs (i.e. pk free ones)
	// This is done similar to s_min_folding::fold_sequence_restricted()
	for (j=0; j < nb_nucleotides; j++)
    {
        for (i=0; i<j; i++)
        {
            // V(i,j) = infinity if i restricted or j restricted and pair of i is not j
            if ((fres[i].pair > -1 && fres[i].pair !=j) || (fres[j].pair > -1 && fres[j].pair != i))
                continue;
            if (fres[i].pair == -1 || fres[j].pair == -1)   // i or j MUST be unpaired
                continue;
            V->compute_energy_restricted_pkonly (i, j, fres); // in s_energy_matrix in simfold package


        }
        // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
        VM->compute_energy_WM_restricted_pkonly (j, fres); // added April 18, 2012

    }


	for (j=0; j < nb_nucleotides; j++)
    {
		// Hosna, March 19, 2012
        for (i =j; i >= 0; i--)//for (i=0; i<=j; i++)
        {
			WMB->compute_energies(i,j);

        	vm->WM_compute_energy(i,j);
        }

	}


	// end of addition at March 8, 2012, Hosna

	for (j= 1; j < nb_nucleotides; j++)
    {
    	this->compute_W_restricted_pkonly(j,fres);
    }
    energy = this->W[nb_nucleotides-1]/100.0;
//    printf("energy = %f \n", energy);


    // backtrack
    // first add (0,n-1) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 0;
    stack_interval->j = nb_nucleotides - 1;
    stack_interval->energy = W[nb_nucleotides-1];
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

    while ( cur_interval != NULL)
    {
        stack_interval = stack_interval->next;
        backtrack_restricted_pkonly (cur_interval, fres); // added April 30, 2012 Hosna
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }


    if (debug)
    {
        print_result ();
    }

    destruct_str_features(nb_nucleotides, fres);
    delete [] h_fres;
    delete [] fres;
    return energy;

}


double W_final::hfold(){

	double energy;
    int i, j;

	h_str_features *h_fres;
    if ((h_fres = new h_str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "h_str_features");
    // detect the structure features
    detect_h_structure_features (restricted, h_fres);
    WMB->set_features(h_fres);
    WMB->initialize();

    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    detect_structure_features (restricted, fres);



    // Hosna: June 28, 2007
    // set the features for checking
    v->set_features(fres);

    // Hosna: July 2nd, 2007
    // set the VM matrix for VM_final
    vm->set_VM_matrix(VM);


	// TODO:
	// I think I shoud fill simfold tables here, before filling the HFold tables
	// Hosna, March 8, 2012

	// 1) fill all th ematrices simfold needs (i.e. pk free ones)
	// This is done similar to s_min_folding::fold_sequence_restricted()
	for (j=0; j < nb_nucleotides; j++)
    {
        for (i=0; i<j; i++)
        {
            // V(i,j) = infinity if i restricted or j restricted and pair of i is not j
            if ((fres[i].pair > -1 && fres[i].pair !=j) || (fres[j].pair > -1 && fres[j].pair != i))
                continue;
            if (fres[i].pair == -1 || fres[j].pair == -1)   // i or j MUST be unpaired
                continue;
            V->compute_energy_restricted (i, j, fres);

        }
        // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
        VM->compute_energy_WM_restricted (j, fres);

		// test V values
		/*
		 for (i=0; i<j; i++)
		 {
		 if (fres[i].pair ==j && fres[j].pair ==i){
		 printf("---->> V(%d,%d) = %d \n",i,j, V->get_energy(i,j));

		 }
		 }
		 */
    }

	for (j=0; j < nb_nucleotides; j++)
    {
        for (i =j; i >= 0; i--)//for (i=0; i<=j; i++)
        {
			WMB->compute_energies(i,j);

			vm->WM_compute_energy(i,j);
			//        	if (debug){
			//        		printf("WM_final(%d,%d) = %d \n",i,j,vm->get_energy_WM(i,j));
			//        	}
        }

	}


	// end of addition at March 8, 2012, Hosna

	for (j= 1; j < nb_nucleotides; j++)
    {
    	this->compute_W_restricted(j,fres);
    }
    energy = this->W[nb_nucleotides-1]/100.0;
	//    printf("energy = %f \n", energy);





    // backtrack
    // first add (0,n-1) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 0;
    stack_interval->j = nb_nucleotides - 1;
    stack_interval->energy = W[nb_nucleotides-1];
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

    while ( cur_interval != NULL)
    {
        stack_interval = stack_interval->next;
        backtrack_restricted (cur_interval, fres);
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }


    if (debug)
    {
        print_result ();
    }

    destruct_str_features(nb_nucleotides, fres);
    delete [] h_fres;
    delete [] fres;
    return energy;

}


//AP
double W_final::hfold_emodel() { //kevin debug
	int KEVIN_DEBUG = 0;

	double energy;
    int i, j;

	h_str_features *h_fres;
    if ((h_fres = new h_str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "h_str_features");
    // detect the structure features
    detect_h_structure_features (restricted, h_fres);
    WMB->set_features(h_fres);
    WMB->initialize();

    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    detect_structure_features (restricted, fres);

    // Hosna: June 28, 2007
    // set the features for checking
    v->set_features(fres);

    // Hosna: July 2nd, 2007
    // set the VM matrix for VM_final
    vm->set_VM_matrix(VM);

	// TODO:
	// I think I shoud fill simfold tables here, before filling the HFold tables
	// Hosna, March 8, 2012

	// 1) fill all th ematrices simfold needs (i.e. pk free ones)
	// This is done similar to s_min_folding::fold_sequence_restricted()
	for (j=1; j < nb_nucleotides; j++)
    {
        for (i=0; i<j; i++)
        {
            // V(i,j) = infinity if i restricted or j restricted and pair of i is not j
            if ((fres[i].pair > -1 && fres[i].pair !=j) || (fres[j].pair > -1 && fres[j].pair != i))
                continue;
            if (fres[i].pair == -1 || fres[j].pair == -1)   // i or j MUST be unpaired
                continue;
			//16 Aug 2017 kevin
			//added fourth argument to check if ij is weakly closed
            V->compute_energy_restricted_emodel (i, j, fres,WMB->is_weakly_closed(i,j));

        }
		
        // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
        VM->compute_energy_WM_restricted_emodel (j, fres, energy_models);
	
                 
		// test V values
		/*
		 for (i=0; i<j; i++)
		 {
		 if (fres[i].pair ==j && fres[j].pair ==i){
		 printf("---->> V(%d,%d) = %d \n",i,j, V->get_energy(i,j));

		 }
		 }
		 */
    }



	for (j=1; j < nb_nucleotides; j++) {
        for (i =j; i >= 0; i--) {//for (i=0; i<=j; i++) {
			    WMB->compute_energies_emodel(i,j,energy_models); // TODO need to check this one
			    vm->WM_compute_energy(i,j); 
        }
	}

//exit(999);

	// end of addition at March 8, 2012, Hosna
	for (j= 1; j < nb_nucleotides; j++) {
    	this->compute_W_restricted_emodel(j,fres);
	}
	if(KEVIN_DEBUG){
		//printf("fres:\n");
		//for (i=0; i < nb_nucleotides; i++){printf("%c i=%d f_type=%c f_pair=%d\n",restricted[i],i,fres[i].type, fres[i].pair);} //kevin debug
		printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~W1:\n");
		for (i=0; i < nb_nucleotides; i++){printf("%c i=%d W=%d\n",restricted[i],i,W[i]);} //kevin debug
	}

    // backtrack
    // first add (0,n-1) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 0;
    stack_interval->j = nb_nucleotides - 1;
    stack_interval->energy = W[nb_nucleotides-1];
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

/*
	//kevin debug
	if(KEVIN_DEBUG){
		printf("before backtrack_restricted_emodel\nstructure: %s\n",structure);
    	for (i=0; i < nb_nucleotides; i++){printf("%c i=%d f_type=%c f_pair=%d\n",restricted[i],i,fres[i].type, fres[i].pair);} //kevin debug
	}
*/

//printf("before backtrack\n");
    while ( cur_interval != NULL) {
        stack_interval = stack_interval->next;
        backtrack_restricted_emodel (cur_interval, fres); // TODO do we need to check this one?
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }
//printf("end of backtrack\n");
	// The energy calculation is now placed after backtrack is run because we need the contents of f[] (aka typedef struct minimum_fold) in order to determine if the final structure is pseudoknoted or not. If it is then we add the start_hybrid_penalty to our final energy and divide it by 100.
    energy = this->W[nb_nucleotides-1]; //nb_nucleotides-1

	// Ian Wark and Kevin July 20 2017
	// We don't need this anymore. It is now in s_energy_matrix::compute_energy_restricted_emodel after calculating hairpin
	/*
	for (i = 0; i < linker_pos; i++) {
		if (f[i].pair > linker_pos+linker_length-1) {
			energy += start_hybrid_penalty;
			break;
		}
	}
	*/

	energy /= 100.0;

	if (debug)
    {
        print_result ();
    }

    destruct_str_features(nb_nucleotides, fres);
    delete [] h_fres;
    delete [] fres;

    return energy;
}

//kevin 18 July
double W_final::call_simfold_emodel(){
	double energy;
    int i, j;

	//16 Aug 2017 kevin
	//added this block so we can use WMB->weakly closed in V->compute_energy_restricted_emodel 
	h_str_features *h_fres;
    if ((h_fres = new h_str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "h_str_features");
    // detect the structure features
    detect_h_structure_features (restricted, h_fres);
    WMB->set_features(h_fres);
    WMB->initialize();

    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    detect_structure_features (restricted, fres);

    /*
    for (i=0; i < nb_nucleotides; i++)
        if (fres[i].pair != -1)
            printf ("%d pairs %d, type %c\n", i, fres[i].pair, fres[i].type);
    */

    for (j=1; j < nb_nucleotides; j++)
    {
        for (i=0; i<j; i++)
        {
            // V(i,j) = infinity if i restricted or j restricted and pair of i is not j
            if ((fres[i].pair > -1 && fres[i].pair !=j) || (fres[j].pair > -1 && fres[j].pair != i))
                continue;
            if (fres[i].pair == -1 || fres[j].pair == -1)   // i or j MUST be unpaired
                continue;
			//16 Aug 2017 kevin
			//added fourth argument to check if ij is weakly closed
            V->compute_energy_restricted_emodel (i, j, fres,WMB->is_weakly_closed(i,j));
        }
        // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
        VM->compute_energy_WM_restricted_emodel (j, fres, energy_models);
    }
    for (j=1; j < nb_nucleotides; j++)
    {
        compute_W_restricted_simfold_emodel (j, fres);
    }
    energy = W[nb_nucleotides-1]/100.0;

    if (debug)
    {
        for (j=1; j < nb_nucleotides; j++)
        {
            printf ("W(%d) = %d\n", j, W[j]);
        }
    }

    // backtrack
    // first add (0,n-1) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 0;
    stack_interval->j = nb_nucleotides - 1;
    stack_interval->energy = W[nb_nucleotides-1];
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

    while ( cur_interval != NULL)
    {
        stack_interval = stack_interval->next;
        backtrack_restricted_simfold_emodel (cur_interval, fres);
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }

    if (debug)
    {
        print_result ();
    }
    destruct_str_features(nb_nucleotides, fres);
    delete [] fres;
    //delete stack_interval;
    return energy;
}

/*
void W_final::call_simfold_emodel(){
	//kevin todo change theese to be actually multi model
	char config_file[200];
    strcpy (config_file, "./simfold/params/multirnafold.conf");
    double temperature;
    temperature = 37;
	int dna_or_rna;
	dna_or_rna = RNA;
	init_data ("./simfold", config_file, dna_or_rna, temperature);
	fill_data_structures_with_new_parameters ("./simfold/params/turner_parameters_fm363_constrdangles.txt");
	fill_data_structures_with_new_parameters ("./simfold/params/parameters_DP09_chopped.txt");
	double energy = simfold_restricted (sequence, restricted, structure);
	//printf("%s %s %s\n",sequence, restricted, structure);
}
*/


//AP
double W_final::hfold_pkonly_emodel(){
	double energy;
    int i, j;

	h_str_features *h_fres;
    if ((h_fres = new h_str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "h_str_features");
    // detect the structure features
    detect_h_structure_features (restricted, h_fres);
    WMB->set_features(h_fres);
    WMB->initialize();

    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    detect_structure_features (restricted, fres);

    // Hosna: June 28, 2007
    // set the features for checking
    v->set_features(fres);

    // Hosna: July 2nd, 2007
    // set the VM matrix for VM_final
    vm->set_VM_matrix(VM);

	// TODO:
	// I think I shoud fill simfold tables here, before filling the HFold tables
	// Hosna, March 8, 2012

	// 1) fill all th ematrices simfold needs (i.e. pk free ones)
	// This is done similar to s_min_folding::fold_sequence_restricted()
	for (j=1; j < nb_nucleotides; j++) {
        for (i=0; i<j; i++) {
            // V(i,j) = infinity if i restricted or j restricted and pair of i is not j
            if ((fres[i].pair > -1 && fres[i].pair !=j) || (fres[j].pair > -1 && fres[j].pair != i))
                continue;
            if (fres[i].pair == -1 || fres[j].pair == -1)   // i or j MUST be unpaired
                continue;
			//16 Aug 2017 kevin
			//added fourth argument to check if ij is weakly closed
            V->compute_energy_restricted_pkonly_emodel (i, j, fres,WMB->is_weakly_closed(i,j));
        }

	    // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
	    VM->compute_energy_WM_restricted_pkonly_emodel (j, fres, energy_models);
    }

	for (j=1; j < nb_nucleotides; j++) {
		// Hosna, March 19, 2012
        for (i =j; i >= 0; i--) {
			//
			WMB->compute_energies_emodel(i,j,energy_models);

        	vm->WM_compute_energy(i,j);
        }
	}

	// end of addition at March 8, 2012, Hosna

	for (j= 1; j < nb_nucleotides; j++)
    {
    	this->compute_W_restricted_pkonly_emodel(j,fres);
    }

    // backtrack
    // first add (0,n-1) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 0;
    stack_interval->j = nb_nucleotides - 1;
    stack_interval->energy = W[nb_nucleotides-1];
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

	// energy = this->W[nb_nucleotides-1];

    while ( cur_interval != NULL)
    {
        stack_interval = stack_interval->next;
        backtrack_restricted_pkonly_emodel (cur_interval, fres);//TODO
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }


	// The energy calculation is now placed after backtrack is run because we need the contents of f[] (aka typedef struct minimum_fold) in order to determine if the final structure is pseudoknoted or not. If it is then we add the start_hybrid_penalty to our final energy and divide it by 100.
	energy = this->W[nb_nucleotides-1];

	// Ian Wark and Kevin July 20 2017
	// We don't need this anymore. It is now in s_energy_matrix::compute_energy_restricted_pkonly_emodel after calculating hairpin
/*
	for (i = 0; i < linker_pos; i++) {
		if (f[i].pair > linker_pos+linker_length-1) {
			energy += start_hybrid_penalty;
			break;
		}
	}
*/
  	energy /= 100.0;

    destruct_str_features(nb_nucleotides, fres);
    delete [] h_fres;
    delete [] fres;
    return energy;
}

void W_final::return_structure(char *structure){
	strcpy (structure, this->structure);
	//s_min_folding::return_structure(structure);
}

void W_final::compute_W_restricted (int j, str_features *fres)
// compute W(j)
{
    int m1, m2, m3;
    int must_choose_this_branch;
    m1 = W[j-1];
    m2 = compute_W_br2_restricted (j, fres, must_choose_this_branch);
    m3 = compute_W_br3_restricted (j, fres);
    if (WMB->is_weakly_closed(0,j) < 0){
    	W[j] = INF;
    	return;
    }

    if (must_choose_this_branch)
    {
        W[j] = MIN(m2,m3);
    }
    else
    {
        W[j] = MIN(m1,MIN(m2,m3));
    }
}

//AP
void W_final::compute_W_restricted_emodel (int j, str_features *fres)
// compute W(j)
{
    int m1 = W[j-1];
	int m2, m3;
    int must_choose_this_branch = 0;

	//AP. This was moved to be checked before the m2 and m3 calculation to save time. Previously it was after the calculation.
	if (WMB->is_weakly_closed(0,j) < 0){
		W[j] = INF;
		return;
	}

    m2 = compute_W_br2_restricted_emodel (j, fres, must_choose_this_branch);
    m3 = compute_W_br3_restricted (j, fres); // TODO ian does this need an emodel?

	if (must_choose_this_branch) {
		W[j] = MIN(m2,m3);
	} else {
		W[j] = MIN(m1,MIN(m2,m3));
	}
}

//AP
//Kevin July 24 2017
void W_final::compute_W_restricted_simfold_emodel (int j, str_features *fres)
// compute W(j)
{
    int m1 = W[j-1];
	int m2;
    int must_choose_this_branch = 0;

    m2 = compute_W_br2_restricted_simfold_emodel (j, fres, must_choose_this_branch);
	if (must_choose_this_branch) {
		W[j] = m2;
	} else {
		W[j] = MIN(m1,m2);
	}
}

//AP
void W_final::compute_W_restricted_pkonly_emodel (int j, str_features *fres)
// compute W(j)
{
    int m1 = W[j-1];
	int m2, m3;
    int must_choose_this_branch = 0;

	//AP. This was moved to be checked before the m2 and m3 calculation to save time. Previously it was after the calculation.
	if (WMB->is_weakly_closed(0,j) < 0){
		W[j] = INF;
		return;
	}

    m2 = compute_W_br2_restricted_pkonly_emodel (j, fres, must_choose_this_branch);
    m3 = compute_W_br3_restricted (j, fres);

	if (must_choose_this_branch) {
		W[j] = MIN(m2,m3);
	} else {
		W[j] = MIN(m1,MIN(m2,m3));
	}
}



void W_final::compute_W_restricted_pkonly (int j, str_features *fres)
// compute W(j)
{
    int m1, m2, m3;
    int must_choose_this_branch;
    m1 = W[j-1];
    m2 = compute_W_br2_restricted_pkonly (j, fres, must_choose_this_branch);
    m3 = compute_W_br3_restricted (j, fres);
    if (WMB->is_weakly_closed(0,j) < 0){
    	W[j] = INF;
    	return;
    }

    if (must_choose_this_branch)
    {
        W[j] = MIN(m2,m3);
    }
    else
    {
        W[j] = MIN(m1,MIN(m2,m3));
    }

}


int W_final::compute_W_br2_restricted (int j, str_features *fres, int &must_choose_this_branch)
{
	int min = INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;

	must_choose_this_branch = 0;
    for (i=0; i<=j-1; i++)    // TURN shouldn't be there
    {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

        // We don't need to make sure i and j don't have to pair with something else,
        //  because that would be INF - done in fold_sequence_restricted
        acc = (i-1>0) ? W[i-1]: 0;

        energy_ij = v->get_energy(i,j);

        if (energy_ij < INF)
        {
            tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
            if (tmp < min)
            {
                min = tmp;
                chosen = 21;        best_i = i;
                if (fres[i].pair == j){
					must_choose_this_branch = 1;
				}
                else                    must_choose_this_branch = 0;
            }
        }

        // I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
        if (fres[i].pair <= -1 && i+1 < j)
        {
            energy_ij = v->get_energy(i+1,j);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
                PARAMTYPE dan = dangle_bot [int_sequence[j]]
											[int_sequence[i+1]]
											[int_sequence[i]];
				//Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i+1].pair == j && fres[j].pair==i+1){
					dan = MIN(0,dan);
				}
				tmp += dan;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 22;  best_i = i;
                    if (fres[i+1].pair == j){
						must_choose_this_branch = 1;
					}
                    else                      must_choose_this_branch = 0;
                }
            }
        }

        // I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
        if (fres[j].pair <= -1 && i < j-1)
        {
            energy_ij = v->get_energy(i,j-1);
            if (energy_ij < INF)
            {
				PARAMTYPE AU_pen=AU_penalty (int_sequence[i],int_sequence[j-1]);
                tmp = energy_ij + AU_pen+ acc;
				PARAMTYPE dan = dangle_top  [int_sequence [j-1]]
											[int_sequence [i]]
											[int_sequence [j]];
				//Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i].pair == j-1 && fres[j-1].pair==i){
					dan = MIN(0,dan);
				 }
                tmp += dan;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 23;  best_i = i;
                    if (fres[i].pair == j-1){
						must_choose_this_branch = 1;
					}
                    else                      must_choose_this_branch = 0;
                }
            }
        }

        if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
        {
            energy_ij = v->get_energy(i+1,j-1);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
				PARAMTYPE dan_bot = dangle_bot [int_sequence[j-1]]
												[int_sequence[i+1]]
												[int_sequence[i]];

				PARAMTYPE dan_top = dangle_top [int_sequence [j-1]]
									[int_sequence [i+1]]
									[int_sequence [j]];
                //Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i+1].pair == j-1 && fres[j-1].pair==i+1){
					dan_bot = MIN(0,dan_bot);
					dan_top = MIN(0,dan_top);
				}
				tmp += dan_bot;
                tmp += dan_top;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 24;  best_i = i;
                    if (fres[i+1].pair == j-1){
						must_choose_this_branch = 1;
					}
                    else                        must_choose_this_branch = 0;
                }
            }
        }
    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min;
}

//kevin
int W_final::compute_W_br2_restricted_simfold_emodel (int j, str_features *fres, int &must_choose_this_branch) {
	int min = INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;
	energy_model *model;

	// Ian Wark and Kevin July 20 2017
	// If j or j-1 is X it cannot be paired
	// j is done here to save time if it is invalid
	if (int_sequence[j] == X || int_sequence[j-1] == X)
        return INF;


	must_choose_this_branch = 0;
    for (i=0; i<=j-1; i++) {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

		//AP
		//if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
		//	continue;

		// Ian Wark and Kevin July 20 2017
        // If i or i+1 is X it cannot be paired
        // i is done here because it needs to be in the for
        if (int_sequence[i] == X || int_sequence[i+1] == X )
			continue;

		for (auto &energy_model : *energy_models) {
			model = &energy_model;
			model->energy_value = min;
		    // We don't need to make sure i and j don't have to pair with something else,
		    //  because that would be INF - done in fold_sequence_restricted
		    acc = (i-1>0) ? W[i-1]: 0;

		    energy_ij = V->get_energy(i,j);
/*
			if (energy_ij != INF && j == 23) {//if (energy_ij != INF && j == 165) {
				acc = (i-1>0) ? W[i-1]: 0;
			}
*/
		    if (energy_ij < INF)
		    {
		        tmp = energy_ij + AU_penalty_emodel (int_sequence[i], int_sequence[j], model) + acc;
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		            chosen = 21;
					best_i = i;
		            if (fres[i].pair == j){
						must_choose_this_branch = 1;
					} else {
						must_choose_this_branch = 0;
					}
				}
		    }

		    // I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
		    if (fres[i].pair <= -1 && i+1 < j)
		    {
		        energy_ij = V->get_energy(i+1,j);
		        if (energy_ij < INF)
		        {
		            tmp = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j], model) + acc;
		            PARAMTYPE dan = model->dangle_bot [int_sequence[j]]
												[int_sequence[i+1]]
												[int_sequence[i]];
					//Hosna, March 27, 2012
					// dangle is INF if the bases are non-canonical and \leq 0 otherwise
					// to accommodate non-canonical base pairing in input structure I add the MIN
					if (fres[i+1].pair == j && fres[j].pair==i+1){
						dan = MIN(0,dan);
					}
					tmp += dan;
		            if (tmp < model->energy_value)
		            {
		                model->energy_value = tmp;
		                chosen = 22;
						best_i = i;
		                if (fres[i+1].pair == j){
							must_choose_this_branch = 1;
						} else {
							must_choose_this_branch = 0;
		            	}
					}
		        }
		    }

		    // I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
		    if (fres[j].pair <= -1 && i < j-1)
		    {
		        energy_ij = V->get_energy(i,j-1);
		        if (energy_ij < INF)
		        {
					PARAMTYPE AU_pen = AU_penalty_emodel (int_sequence[i], int_sequence[j-1], model);
		            tmp = energy_ij + AU_pen+ acc;
					PARAMTYPE dan = model->dangle_top  [int_sequence [j-1]]
												[int_sequence [i]]
												[int_sequence [j]];
					//Hosna, March 27, 2012
					// dangle is INF if the bases are non-canonical and \leq 0 otherwise
					// to accommodate non-canonical base pairing in input structure I add the MIN
					if (fres[i].pair == j-1 && fres[j-1].pair==i){
						dan = MIN(0,dan);
					 }
		            tmp += dan;
		            if (tmp < model->energy_value)
		            {
		                model->energy_value = tmp;
		                chosen = 23;
						best_i = i;
		                if (fres[i].pair == j-1){
							must_choose_this_branch = 1;
						} else {
							must_choose_this_branch = 0;
						}
					}
		        }
		    }

		    if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
		    {
		        energy_ij = V->get_energy(i+1,j-1);
		        if (energy_ij < INF)
		        {
		            tmp = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j-1], model) + acc;
					PARAMTYPE dan_bot = model->dangle_bot [int_sequence[j-1]]
													[int_sequence[i+1]]
													[int_sequence[i]];

					PARAMTYPE dan_top = model->dangle_top [int_sequence [j-1]]
										[int_sequence [i+1]]
										[int_sequence [j]];
		            //Hosna, March 27, 2012
					// dangle is INF if the bases are non-canonical and \leq 0 otherwise
					// to accommodate non-canonical base pairing in input structure I add the MIN
					if (fres[i+1].pair == j-1 && fres[j-1].pair==i+1){
						dan_bot = MIN(0,dan_bot);
						dan_top = MIN(0,dan_top);
					}
					tmp += dan_bot;
		            tmp += dan_top;
		            if (tmp < model->energy_value)
		            {
		                model->energy_value = tmp;
		                chosen = 24;
						best_i = i;
		                if (fres[i+1].pair == j-1){
							must_choose_this_branch = 1;
						} else {
							must_choose_this_branch = 0;
		            	}
					}
		        }
		    }
		}
		tmp = emodel_energy_function (i, j, energy_models);
		if (tmp != INF) {
			min = tmp;
		}

    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min;
}

//AP
int W_final::compute_W_br2_restricted_emodel (int j, str_features *fres, int &must_choose_this_branch) {
	int min = INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;
	energy_model *model;

	// Ian Wark and Kevin July 20 2017
	// If j or j-1 is X it cannot be paired
	// j is done here to save time if it is invalid
	if (int_sequence[j] == X || int_sequence[j-1] == X)
        return INF;


	must_choose_this_branch = 0;
    for (i=0; i<=j-1; i++) {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

		//AP
		//if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
		//	continue;

		// Ian Wark and Kevin July 20 2017
        // If i or i+1 is X it cannot be paired
        // i is done here because it needs to be in the for
        if (int_sequence[i] == X || int_sequence[i+1] == X )
			continue;

		for (auto &energy_model : *energy_models) {
			model = &energy_model;
			model->energy_value = min;
		    // We don't need to make sure i and j don't have to pair with something else,
		    //  because that would be INF - done in fold_sequence_restricted
		    acc = (i-1>0) ? W[i-1]: 0;

		    energy_ij = v->get_energy(i,j);
			
/*
			if (energy_ij != INF && j == 23) {//if (energy_ij != INF && j == 165) {
				acc = (i-1>0) ? W[i-1]: 0;
			}
*/
		    if (energy_ij < INF)
		    {
		        tmp = energy_ij + AU_penalty_emodel (int_sequence[i], int_sequence[j], model) + acc;
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		            chosen = 21;
					best_i = i;
		            if (fres[i].pair == j){
						must_choose_this_branch = 1;
					} else {
						must_choose_this_branch = 0;
					}
				}
		    }

		    // I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
		    if (fres[i].pair <= -1 && i+1 < j)
		    {
		        energy_ij = v->get_energy(i+1,j);
		        if (energy_ij < INF)
		        {
		            tmp = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j], model) + acc;
		            PARAMTYPE dan = model->dangle_bot [int_sequence[j]]
												[int_sequence[i+1]]
												[int_sequence[i]];
					//Hosna, March 27, 2012
					// dangle is INF if the bases are non-canonical and \leq 0 otherwise
					// to accommodate non-canonical base pairing in input structure I add the MIN
					if (fres[i+1].pair == j && fres[j].pair==i+1){
						dan = MIN(0,dan);
					}
					tmp += dan;
		            if (tmp < model->energy_value)
		            {
		                model->energy_value = tmp;
		                chosen = 22;
						best_i = i;
		                if (fres[i+1].pair == j){
							must_choose_this_branch = 1;
						} else {
							must_choose_this_branch = 0;
		            	}
					}
		        }
		    }

		    // I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
		    if (fres[j].pair <= -1 && i < j-1)
		    {
		        energy_ij = v->get_energy(i,j-1);
		        if (energy_ij < INF)
		        {
					PARAMTYPE AU_pen = AU_penalty_emodel (int_sequence[i], int_sequence[j-1], model);
		            tmp = energy_ij + AU_pen+ acc;
					PARAMTYPE dan = model->dangle_top  [int_sequence [j-1]]
												[int_sequence [i]]
												[int_sequence [j]];
					//Hosna, March 27, 2012
					// dangle is INF if the bases are non-canonical and \leq 0 otherwise
					// to accommodate non-canonical base pairing in input structure I add the MIN
					if (fres[i].pair == j-1 && fres[j-1].pair==i){
						dan = MIN(0,dan);
					 }
		            tmp += dan;
		            if (tmp < model->energy_value)
		            {
		                model->energy_value = tmp;
		                chosen = 23;
						best_i = i;
		                if (fres[i].pair == j-1){
							must_choose_this_branch = 1;
						} else {
							must_choose_this_branch = 0;
						}
					}
		        }
		    }

		    if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
		    {
		        energy_ij = v->get_energy(i+1,j-1);
		        if (energy_ij < INF)
		        {
		            tmp = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j-1], model) + acc;
					PARAMTYPE dan_bot = model->dangle_bot [int_sequence[j-1]]
													[int_sequence[i+1]]
													[int_sequence[i]];

					PARAMTYPE dan_top = model->dangle_top [int_sequence [j-1]]
										[int_sequence [i+1]]
										[int_sequence [j]];
		            //Hosna, March 27, 2012
					// dangle is INF if the bases are non-canonical and \leq 0 otherwise
					// to accommodate non-canonical base pairing in input structure I add the MIN
					if (fres[i+1].pair == j-1 && fres[j-1].pair==i+1){
						dan_bot = MIN(0,dan_bot);
						dan_top = MIN(0,dan_top);
					}
					tmp += dan_bot;
		            tmp += dan_top;
		            if (tmp < model->energy_value)
		            {
		                model->energy_value = tmp;
		                chosen = 24;
						best_i = i;
		                if (fres[i+1].pair == j-1){
							must_choose_this_branch = 1;
						} else {
							must_choose_this_branch = 0;
		            	}
					}
		        }
		    }
		}
		tmp = emodel_energy_function (i, j, energy_models);
		if (tmp != INF) {
			min = tmp;
		}
		
    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min;
}

//AP
int W_final::compute_W_br2_restricted_pkonly_emodel (int j, str_features *fres, int &must_choose_this_branch) {
	int min = INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;
	energy_model *model;

  // Ian Wark and Kevin July 20 2017
	// If j or j-1 is X it cannot be paired
	// j is done here to save time if it is invalid
	if (int_sequence[j] == X || int_sequence[j-1] == X)
        return INF;

	must_choose_this_branch = 0;
    for (i=0; i<=j-1; i++) {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

		//AP
		//if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
		//	continue;

		// Ian Wark and Kevin July 20 2017
    // If i or i+1 is X it cannot be paired
    // i is done here because it needs to be in the for
    if (int_sequence[i] == X || int_sequence[i+1] == X )
			continue;

		for (auto &energy_model : *energy_models) {
			model = &energy_model;
			model->energy_value = min;
		    // We don't need to make sure i and j don't have to pair with something else,
		    //  because that would be INF - done in fold_sequence_restricted
		    acc = (i-1>0) ? W[i-1]: 0;

		    energy_ij = (fres[j].pair == i && fres[i].pair ==j)? v->get_energy(i,j) : INF; //v->get_energy(i,j);

		    if (energy_ij < INF)
		    {
		        tmp = energy_ij + AU_penalty_emodel (int_sequence[i], int_sequence[j], model) + acc;
		        if (tmp < model->energy_value)
		        {
		            model->energy_value = tmp;
		            chosen = 21;
					best_i = i;
		            if (fres[i].pair == j){
						must_choose_this_branch = 1;
					} else {
						must_choose_this_branch = 0;
					}
		        }
		    }

		    // I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
			// in the pk_only version we don't add any pseudoknot free base pairs
		    if (fres[i].pair <= -1 && i+1 < j && fres[i+1].pair ==j && fres[j].pair==i+1)
		    {
		        energy_ij = v->get_energy(i+1,j);
		        if (energy_ij < INF)
		        {
		            tmp = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j], model) + acc;
		            PARAMTYPE dan = model->dangle_bot [int_sequence[j]]
													[int_sequence[i+1]]
													[int_sequence[i]];
					//Hosna, March 27, 2012
					// dangle is INF if the bases are non-canonical and \leq 0 otherwise
					// to accommodate non-canonical base pairing in input structure I add the MIN
					if (fres[i+1].pair == j && fres[j].pair==i+1){
						dan = MIN(0,dan);
					}
					tmp += dan;
		            if (tmp < model->energy_value)
		            {
		                model->energy_value = tmp;
		                chosen = 22;
						best_i = i;
		                if (fres[i+1].pair == j){
							must_choose_this_branch = 1;
						} else {
							must_choose_this_branch = 0;
						}
		            }
		        }
		    }

		    // I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
			// in the pkonly version we don't add any pseudoknot free base pairs
		    if (fres[j].pair <= -1 && i < j-1 && fres[i].pair ==j-1 && fres[j-1].pair ==i)
		    {
		        energy_ij = v->get_energy(i,j-1);
		        if (energy_ij < INF)
		        {
					PARAMTYPE AU_pen = AU_penalty_emodel (int_sequence[i], int_sequence[j-1], model);
		            tmp = energy_ij + AU_pen+ acc;
					PARAMTYPE dan = model->dangle_top [int_sequence [j-1]]
													[int_sequence [i]]
													[int_sequence [j]];
					//Hosna, March 27, 2012
					// dangle is INF if the bases are non-canonical and \leq 0 otherwise
					// to accommodate non-canonical base pairing in input structure I add the MIN
					if (fres[i].pair == j-1 && fres[j-1].pair==i){
						dan = MIN(0,dan);
					}
		            tmp += dan;
		            if (tmp < model->energy_value)
		            {
		                model->energy_value = tmp;
		                chosen = 23;
						best_i = i;
		                if (fres[i].pair == j-1){
							must_choose_this_branch = 1;
						} else {
							must_choose_this_branch = 0;
						}
		            }
		        }
		    }

			// in the pkonly version we don't add any pseudoknot free base pairs
		    if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1 && fres[i+1].pair==j-1 && fres[j-1].pair==i+1)
		    {
		        energy_ij = v->get_energy(i+1,j-1);
		        if (energy_ij < INF)
		        {
		            tmp = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j-1], model) + acc;
					PARAMTYPE dan_bot = model->dangle_bot [int_sequence[j-1]]
														[int_sequence[i+1]]
														[int_sequence[i]];

					PARAMTYPE dan_top = model->dangle_top [int_sequence [j-1]]
														[int_sequence [i+1]]
														[int_sequence [j]];
		            //Hosna, March 27, 2012
					// dangle is INF if the bases are non-canonical and \leq 0 otherwise
					// to accommodate non-canonical base pairing in input structure I add the MIN
					if (fres[i+1].pair == j-1 && fres[j-1].pair==i+1){
						dan_bot = MIN(0,dan_bot);
						dan_top = MIN(0,dan_top);
					}
					tmp += dan_bot;
		            tmp += dan_top;
		            if (tmp < model->energy_value)
		            {
		                model->energy_value = tmp;
		                chosen = 24;
						best_i = i;
		                if (fres[i+1].pair == j-1){
							must_choose_this_branch = 1;
						} else {
							must_choose_this_branch = 0;
						}
		            }
		        }
		    }
		}
		tmp = emodel_energy_function (i, j, energy_models);
		if (tmp != INF) {
			min = tmp;
		}

    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min;

}

int W_final::compute_W_br2_restricted_pkonly (int j, str_features *fres, int &must_choose_this_branch)
{
	int min = INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;

	must_choose_this_branch = 0;
    for (i=0; i<=j-1; i++)    // TURN shouldn't be there
    {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

        // We don't need to make sure i and j don't have to pair with something else,
        //  because that would be INF - done in fold_sequence_restricted
        acc = (i-1>0) ? W[i-1]: 0;

        energy_ij = (fres[j].pair == i && fres[i].pair ==j)? v->get_energy(i,j) : INF; //v->get_energy(i,j);

        if (energy_ij < INF)
        {
            tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
            if (tmp < min)
            {
                min = tmp;
                chosen = 21;        best_i = i;
                if (fres[i].pair == j){
					must_choose_this_branch = 1;
				}
                else                    must_choose_this_branch = 0;
            }
        }

        // I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
		// in the pk_only version we don't add any pseudoknot free base pairs
        if (fres[i].pair <= -1 && i+1 < j && fres[i+1].pair ==j && fres[j].pair==i+1)
        {
            energy_ij = v->get_energy(i+1,j);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
                PARAMTYPE dan = dangle_bot [int_sequence[j]]
				[int_sequence[i+1]]
				[int_sequence[i]];
				//Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i+1].pair == j && fres[j].pair==i+1){
					dan = MIN(0,dan);
				}
				tmp += dan;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 22;  best_i = i;
                    if (fres[i+1].pair == j){
						must_choose_this_branch = 1;
					}
                    else                      must_choose_this_branch = 0;
                }
            }
        }

        // I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
		// in the pkonly version we don't add any pseudoknot free base pairs
        if (fres[j].pair <= -1 && i < j-1 && fres[i].pair ==j-1 && fres[j-1].pair ==i)
        {
            energy_ij = v->get_energy(i,j-1);
            if (energy_ij < INF)
            {
				PARAMTYPE AU_pen=AU_penalty (int_sequence[i],int_sequence[j-1]);
                tmp = energy_ij + AU_pen+ acc;
				PARAMTYPE dan = dangle_top  [int_sequence [j-1]]
				[int_sequence [i]]
				[int_sequence [j]];
				//Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i].pair == j-1 && fres[j-1].pair==i){
					dan = MIN(0,dan);
				}
                tmp += dan;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 23;  best_i = i;
                    if (fres[i].pair == j-1){
						must_choose_this_branch = 1;
					}
                    else                      must_choose_this_branch = 0;
                }
            }
        }

		// in the pkonly version we don't add any pseudoknot free base pairs
        if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1 && fres[i+1].pair==j-1 && fres[j-1].pair==i+1)
        {
            energy_ij = v->get_energy(i+1,j-1);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
				PARAMTYPE dan_bot = dangle_bot [int_sequence[j-1]]
				[int_sequence[i+1]]
				[int_sequence[i]];

				PARAMTYPE dan_top = dangle_top [int_sequence [j-1]]
				[int_sequence [i+1]]
				[int_sequence [j]];
                //Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i+1].pair == j-1 && fres[j-1].pair==i+1){
					dan_bot = MIN(0,dan_bot);
					dan_top = MIN(0,dan_top);
				}
				tmp += dan_bot;
                tmp += dan_top;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 24;  best_i = i;
                    if (fres[i+1].pair == j-1){
						must_choose_this_branch = 1;
					}
                    else                        must_choose_this_branch = 0;
                }
            }
        }
    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min;
}


int W_final::compute_W_br3_restricted(int j, str_features *fres){
	// Hosna June 30, 2007
	// The following would not take care of when
	// we have some unpaired bases before the start of the WMB
	//return WMB->get_energy(0,j) + PS_penalty;
	int min = INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;

//	must_choose_this_branch = 0;
    for (i=0; i<=j-1; i++)    // TURN shouldn't be there
    {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

		// Hosna: July 9, 2007
		// We only chop W to W + WMB when the bases before WMB are free
		if (i == 0 || (WMB->is_weakly_closed(0,i-1) && WMB->is_weakly_closed(i,j))){

	        // We don't need to make sure i and j don't have to pair with something else,
	        //  because that would be INF - done in fold_sequence_restricted
	        acc = (i-1>0) ? W[i-1]: 0;

	        energy_ij = WMB->get_energy(i,j);
	        if (energy_ij < INF)
	        {
	            tmp = energy_ij + PS_penalty + acc;
	            if (tmp < min)
	            {
	                min = tmp;
	                chosen = 31;
	                best_i = i;
	//                if (fres[i].pair == j)  must_choose_this_branch = 1;
	//                else                    must_choose_this_branch = 0;
	            }
	        }

	        // I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
	        if (fres[i].pair <= -1 && i+1 < j)
	        {
	            energy_ij = WMB->get_energy(i+1,j);
	            if (energy_ij < INF)
	            {
	                tmp = energy_ij + PS_penalty + acc;
	//                tmp += dangle_bot [int_sequence[j]]
	//                                [int_sequence[i+1]]
	//                                [int_sequence[i]];
	                if (tmp < min)
	                {
	                    min = tmp;
	                    chosen = 32;
	                    best_i = i;
	//                    if (fres[i+1].pair == j)  must_choose_this_branch = 1;
	//                    else                      must_choose_this_branch = 0;
	                }
	            }
	        }

	        // I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
	        if (fres[j].pair <= -1 && i < j-1)
	        {
	            energy_ij = WMB->get_energy(i,j-1);
	            if (energy_ij < INF)
	            {
	                tmp = energy_ij + PS_penalty + acc;
	//                tmp += dangle_top [int_sequence [j-1]]
	//                                [int_sequence [i]]
	//                                [int_sequence [j]];
	                if (tmp < min)
	                {
	                    min = tmp;
	                    chosen = 33;
	                    best_i = i;
	//                    if (fres[i].pair == j-1)  must_choose_this_branch = 1;
	//                    else                      must_choose_this_branch = 0;
	                }
	            }
	        }

	        if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
	        {
	            energy_ij = WMB->get_energy(i+1,j-1);
	            if (energy_ij < INF)
	            {
	                tmp = energy_ij + PS_penalty + acc;
	//                tmp += dangle_bot [int_sequence[j-1]]
	//                                [int_sequence[i+1]]
	//                                [int_sequence[i]];
	//                tmp += dangle_top [int_sequence [j-1]]
	//                                [int_sequence [i+1]]
	//                                [int_sequence [j]];
	                if (tmp < min)
	                {
	                    min = tmp;
	                    chosen = 34;
	                    best_i = i;
	//                    if (fres[i+1].pair == j-1)  must_choose_this_branch = 1;
	//                    else                        must_choose_this_branch = 0;
	                }
	            }
	        }
		}
    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min;
}

void W_final::backtrack_restricted(seq_interval *cur_interval, str_features *fres){
    char type;

	//Hosna, March 8, 2012
	// changing nested if to switch for optimality
	switch (cur_interval->type){
		case LOOP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i >= j)
				return;
			f[i].pair = j;
			f[j].pair = i;

			// Hosna Jun. 28 2007
			// if the pairing is part of original structure, put '(' and ')' in the structure
			// otherwise make it '[' and ']'
			if (fres[i].pair == j){
				structure[i] = '(';
				structure[j] = ')';
			}else{
				structure[i] = '[';
				structure[j] = ']';
			}

			type = v->get_type (i,j);
			if (debug)
				printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
			// Hosna, March 8, 2012
			// changing nested ifs to switch for optimality
			switch (type){
				case STACK:
			//if (type == STACK)
				{
					f[i].type = STACK;
					f[j].type = STACK;
					if (i+1 < j-1)
						this->insert_node(i+1,j-1, LOOP);
						//insert_node (i+1, j-1, LOOP);
					else
					{
						fprintf (stderr, "NOT GOOD RESTR STACK, i=%d, j=%d\n", i, j);
						exit (0);
					}

				}
					break;
				case HAIRP:
			//else if (type == HAIRP)
				{
					f[i].type = HAIRP;
					f[j].type = HAIRP;
				}
					break;
				case INTER:
			//else if (type == INTER)
				{
					f[i].type = INTER;
					f[j].type = INTER;
					// detect the other closing pair
					int ip, jp, best_ip, best_jp, minq;
					int tmp, min = INF;
					for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1); ip++)
					{
						// Hosna, August 28, 2012
						// TODO: cannot understand why we have th efollowing calculations, as it makes the following case be missed!
						// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
						// (____(_____((__(_______)__))_______________________________)____)
						// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
						// in this example int(5,59,11,27) is falsely missed and is equal to INF
						// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
						//minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
						minq = ip+1;
						for (jp = minq; jp < j; jp++)
						{
							if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
								continue;
							//tmp = VBI->get_energy_str (i,j,ip,jp);
							// Hosna, March 26, 2012
							// modified to accommodate non-canonical base pairing in restricted structure
							tmp = VBI->get_energy_str_restricted(i,j,ip,jp,fres);
							if (tmp < min)
							{
								min = tmp;
								best_ip = ip;
								best_jp = jp;
							}
						}
					}
					if (best_ip < best_jp)
						insert_node (best_ip, best_jp, LOOP);
					else
					{
						fprintf (stderr, "NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
			//else if (type == MULTI)
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					int k, best_k = -1, best_row = -1;
					int tmp= INF, min = INF;
					for (k = i+1; k <= j-1; k++)
					  {
						tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-1);
						if (tmp < min)
						  {
							min = tmp;
							best_k = k;
							best_row = 1;
						  }
						  // TODO:
						  // Hosna, May 1st, 2012
						  // do I need to check for non-canonical base pairings here as well so the dangle values not be INF??
						if (fres[i+1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-1) +
							dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
							misc.multi_free_base_penalty;
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 2;
							}
						}
						if (fres[j-1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-2) +
							dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
							misc.multi_free_base_penalty;
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 3;
							}
						}
						if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-2) +
							dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
							dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
							2*misc.multi_free_base_penalty;
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 4;
							}
						}

						// Hosna: June 28, 2007
						// the last branch of VM, which is WMB_(i+1),(j-1)
						// Hosna: July 5th, 2007
						// add a b_penalty to this case to match the previous cases
						tmp = WMB->get_energy(i+1,j-1)+ a_penalty + PSM_penalty+ b_penalty;
						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
						}

					  }
					switch (best_row)
					  {
					  case 1:
		//              	printf("M_WM(%d,%d) branch 1: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-1);
						insert_node (i+1, best_k, M_WM);
						insert_node (best_k+1, j-1, M_WM);
						break;
					  case 2:
		//              	printf("M_WM(%d,%d) branch 2: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-1);
						insert_node (i+2, best_k, M_WM);
						insert_node (best_k+1, j-1, M_WM);
						break;
					  case 3:
		//              	printf("M_WM(%d,%d) branch 3: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-2);
						insert_node (i+1, best_k, M_WM);
						insert_node (best_k+1, j-2, M_WM);
						break;
					  case 4:
		//              	printf("M_WM(%d,%d) branch 4: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-2);
						insert_node (i+2, best_k, M_WM);
						insert_node (best_k+1, j-2, M_WM);
						break;
					  // Hosna: June 28, 2007
					  // the last branch of VM, which is WMB_(i+1),(j-1)
					  case 5:
		//              	printf("M_WM(%d,%d) branch 3: pushing P_WMB(%d,%d)\n", i,j,i+1,j-1);
						insert_node(i+1,j-1, P_WMB);
						break;
					  }
				}
					break;
			}
		}
			break;
		case FREE:
		{
			int j = cur_interval->j;

			if (j==0) return;

			int min = INF, tmp, best_row, i, best_i, acc, energy_ij;

			if (debug)
				printf ("\t(0,%d) FREE\n", j);

			// this case is for j unpaired, so I have to check that.
			if (fres[j].pair <= -1)
			{
				tmp = W[j-1];

            			if (tmp < min)
            			{
					min = tmp;
					best_row = 0;
				}
			}
			for (i=0; i<=j-1; i++)    // no TURN
			{

				// Don't need to make sure i and j don't have to pair with something else
				//  it's INF, done in fold_sequence_restricted
				acc = (i-1>0) ? W[i-1] : 0;
				energy_ij = v->get_energy(i,j);

				if (energy_ij < INF)
				{
					tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
					if (tmp < min)
					{
					min = tmp;
					best_i = i;
					best_row = 1;
					}
				}

				if (fres[i].pair <= -1)
				{
					energy_ij = v->get_energy(i+1,j);
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
						tmp += dangle_bot [int_sequence[j]]
							[int_sequence[i+1]]
							[int_sequence[i]];
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 2;
						}
					}
				}
				if (fres[j].pair <= -1)
				{
					energy_ij = v->get_energy(i,j-1);
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j-1]) + acc;
						tmp += dangle_top [int_sequence[j-1]]
							[int_sequence[i]]
							[int_sequence[j]];
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 3;
						}
					}
				}
				if (fres[i].pair <= -1 && fres[j].pair <= -1)
				{
					energy_ij = v->get_energy(i+1,j-1);
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
						tmp += dangle_bot [int_sequence[j-1]]
							[int_sequence[i+1]]
							[int_sequence[i]];
						tmp += dangle_top [int_sequence[j-1]]
							[int_sequence[i+1]]
							[int_sequence[j]];
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 4;
						}
					}
				}
			}

			// Hosna: June 28, 2007
			// the last branch of W, which is WMB_i,j
	//        energy_ij = WMB->get_energy(0,j);
	//        if (energy_ij < INF){
	//          	tmp = energy_ij + PS_penalty;
	//           	if (tmp < min){
	//           		min = tmp;
	//           		best_row = 5;
	//           	}
	//        }
		// Hosna June 30, 2007
		// The following would not take care of when
		// we have some unpaired bases before the start of the WMB
		for (i=0; i<=j-1; i++)
		{
			// Hosna: July 9, 2007
			// We only chop W to W + WMB when the bases before WMB are free
			if (i == 0 || (WMB->is_weakly_closed(0,i-1) && WMB->is_weakly_closed(i,j))){

				acc = (i-1>0) ? W[i-1]: 0;

				energy_ij = WMB->get_energy(i,j);

				if (energy_ij < INF)
				{
					tmp = energy_ij + PS_penalty + acc;

					if (tmp < min)
					{
						min = tmp;
						best_row = 5;
						best_i = i;
					}
				}

				// I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
				if (fres[i].pair <= -1 && i+1 < j)
				{
					energy_ij = WMB->get_energy(i+1,j);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 6;
							best_i = i;
						}
					}
				}

				// I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
				if (fres[j].pair <= -1 && i < j-1)
				{
					energy_ij = WMB->get_energy(i,j-1);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 7;
							best_i = i;
						}
					}
				}

				if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
				{
					energy_ij = WMB->get_energy(i+1,j-1);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 8;
							best_i = i;
						}
					}
				}
			}
		}
			switch (best_row)
			{
				case 0:
					//printf("W(%d) case 0: inserting Free (0,%d)\n",j,j-1);
					insert_node (0, j-1, FREE); break;
				case 1:
					//printf("W(%d) case 1: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, LOOP);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 2:
					//printf("W(%d) case 2: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, LOOP);
					if (best_i >= 0)// Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 3:
					//printf("W(%d) case 3: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, LOOP);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 4:
					//printf("W(%d) case 4: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, LOOP);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				// Hosna: June 28, 2007
				// the last branch of W, which is WMB_i,j
				case 5:
					//printf("W(%d) case 5: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, P_WMB);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 6:
					//printf("W(%d) case 6: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 7:
					//printf("W(%d) case 7: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, P_WMB);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 8:
					//printf("W(%d) case 8: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
			}
		}
			break;
		case M_WM:
//  else if(cur_interval->type == M_WM)
		{
			  int i = cur_interval->i;
			  int j = cur_interval->j;
			  int tmp, min = INF;
			  int best_k, best_row;

			  if (debug)
				printf ("\t (%d,%d) M_WM\n", i,j);

			  tmp = v->get_energy(i,j) +
				AU_penalty (int_sequence[i], int_sequence[j]) +
				misc.multi_helix_penalty;
			  if (tmp < min)
				{
				  min = tmp;
				  best_row = 1;
				}
			  if (fres[i].pair <= -1)
			  {
				  tmp = v->get_energy(i+1,j) +
						AU_penalty (int_sequence[i+1], int_sequence[j]) +
						dangle_bot [int_sequence[j]]
						[int_sequence[i+1]]
						[int_sequence[i]] +
						misc.multi_helix_penalty +
						misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 2;
				  }
			  }
			  if (fres[j].pair <= -1)
			  {
				  tmp = v->get_energy(i,j-1) +
						AU_penalty (int_sequence[i], int_sequence[j-1]) +
						dangle_top [int_sequence[j-1]]
									[int_sequence[i]]
									[int_sequence[j]] +
						misc.multi_helix_penalty +
						misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 3;
				  }
			  }
			  if (fres[i].pair <= -1 && fres[j].pair <= -1)
			  {
				  tmp = v->get_energy(i+1,j-1) +
						AU_penalty (int_sequence[i+1], int_sequence[j-1]) +
						dangle_bot [int_sequence[j-1]]
									[int_sequence[i+1]]
									[int_sequence[i]] +
						dangle_top [int_sequence[j-1]]
									[int_sequence[i+1]]
									[int_sequence[j]] +
						misc.multi_helix_penalty +
						2*misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 4;
				  }
			  }
			  if (fres[i].pair <= -1)
			  {
				  tmp = vm->get_energy_WM (i+1,j) + misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 5;
				  }
			  }
			  if (fres[j].pair <= -1)
			  {
				  tmp = vm->get_energy_WM (i,j-1) + misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 6;
				  }
			  }

			  for (int k=i; k < j; k++)
				{
					tmp = vm->get_energy_WM (i, k) + vm->get_energy_WM (k+1, j);
					if (tmp < min)
					  {
						min = tmp;
						best_k = k;
						best_row = 7;
					  }
				}
			  // Hosna: June 28, 2007
			  // the last branch of WW, which is WMB_i,j
			  tmp = WMB->get_energy(i,j)+PSM_penalty;
			  if (tmp < min){
				min = tmp;
				best_row = 8;
			  }

			  switch (best_row)
				{
				  case 1: insert_node (i, j, LOOP); break;
				  case 2: insert_node (i+1, j, LOOP); break;
				  case 3: insert_node (i, j-1, LOOP); break;
				  case 4: insert_node (i+1, j-1, LOOP); break;
				  case 5:
					if (j-i-1 > 0)
					  insert_node (i+1, j, M_WM);
					break;
				  case 6:
					if (j-1-i > 0)
					  insert_node (i, j-1, M_WM);
					break;
				  case 7:
					if (best_k-i > 0)
					  insert_node (i, best_k, M_WM);
					if (j-best_k-1 > 0)
					  insert_node (best_k+1, j, M_WM);
					break;
				  // Hosna: June 28, 2007
				  // the last branch of W, which is WMB_i,j
				  case 8:
					insert_node(i,j,P_WMB);
					break;
				  }
			}
			break;
    // Hosna: Feb 19th 2007
		case P_WMB:
   // else if(cur_interval->type == P_WMB)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WMBP:
    // Hosna: April 18th, 2007
    // changed WMB to case 2 and WMBP
   // else if(cur_interval->type == P_WMBP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VP:
    //else if(cur_interval->type == P_VP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VPP:
    //else if(cur_interval->type == P_VPP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WI:
    //else if(cur_interval->type == P_WI)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_BE:
    //else if(cur_interval->type == P_BE)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WIP:
    //else if(cur_interval->type == P_WIP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		default:
			printf("Should not be here!\n");
	}


}

//AP
void W_final::backtrack_restricted_emodel(seq_interval *cur_interval, str_features *fres){
	int KEVIN_DEBUG = 0;
	if(KEVIN_DEBUG)
		printf("in backtrack_restricted_emodel\n");
    char type;
	energy_model *model;

	//for (auto &energy_model : *energy_models) {
 	//	model = &energy_model;

	//Hosna, March 8, 2012
	// changing nested if to switch for optimality
	switch (cur_interval->type){
		case LOOP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i >= j)
				return;
			f[i].pair = j;
			f[j].pair = i;
		//if(KEVIN_DEBUG){
		//	//AP
		//	if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
		//		return;
		//}
			// Hosna Jun. 28 2007
			// if the pairing is part of original structure, put '(' and ')' in the structure
			// otherwise make it '[' and ']'
			if (fres[i].pair == j){
				structure[i] = '(';
				structure[j] = ')';
			}else{
				structure[i] = '[';
				structure[j] = ']';
			}

			type = v->get_type (i,j);
			if (debug)
				printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
			// Hosna, March 8, 2012
			// changing nested ifs to switch for optimality
			switch (type){
				case STACK:
			//if (type == STACK)
				{
					f[i].type = STACK;
					f[j].type = STACK;
					if (i+1 < j-1)
						this->insert_node(i+1,j-1, LOOP);
						//insert_node (i+1, j-1, LOOP);
					else
					{
						fprintf (stderr, "NOT GOOD RESTR STACK, i=%d, j=%d\n", i, j);
						exit (0);
					}

				}
					break;
				case HAIRP:
			//else if (type == HAIRP)
				{
					f[i].type = HAIRP;
					f[j].type = HAIRP;
				}
					break;
				case INTER:
			//else if (type == INTER)
				{
					f[i].type = INTER;
					f[j].type = INTER;
					// detect the other closing pair
					int ip, jp, best_ip = -1, best_jp = -1, minq;
					int tmp, min = INF;

					//todo kevin confirm
					//23 Aug 2017 kevin and Mahyar
					//variable to store linker_length such that we change the range of ip and jp so we can check if i,ip and jp,j is larger than MAXLOOP (aka 30) properly for the cases where X is between i,j
					//if X is between i,j we treat it as if X does not exist 
					int skip = 0;			
					if(is_cross_model(i,j)){
						skip = linker_length;
					}
					
					
					//todo kevin confirm
					//23 Aug 2017 kevin and Mahyar
					//added +skip to i+MAXLOOP+1
					for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1+skip); ip++)
					{
						// Hosna, August 28, 2012
						// TODO: cannot understand why we have th efollowing calculations, as it makes the following case be missed!
						// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
						// (____(_____((__(_______)__))_______________________________)____)
						// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
						// in this example int(5,59,11,27) is falsely missed and is equal to INF
						// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
						//minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
						
/*
						//todo kevin confirm
						//23 Aug 2017 kevin and Mahyar
						//changed minq to new one (copied from simfold backtrack emodel)
						minq = MAX (j-i+ip-MAXLOOP-2-skip, ip+1);
*/
						minq = ip+1;
						for (jp = minq; jp < j; jp++)
						{
							 
							if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres) ) //|| (fres[ip].pair >= 0 && fres[ip].pair != jp))
								continue;
							
							//tmp = VBI->get_energy_str (i,j,ip,jp);
							// Hosna, March 26, 2012
							// modified to accommodate non-canonical base pairing in restricted structure
							//AP
							for (auto &model : *energy_models) {
		
								//model.energy_value = VBI->get_energy_str_restricted_emodel (i, j, ip, jp, fres, &model);
								//17 Aug 2017 kevin and Mahyar
								//changed the above function call to this one so we dont re-caculate it and just look up the value
								model.energy_value = v->get_energy(ip,jp);
				
							}
							tmp = emodel_energy_function (i, j, energy_models);
					
				
							if (tmp < min)
							{
								min = tmp;
								best_ip = ip;
								best_jp = jp;
							}
						}
					}
					if (best_ip < best_jp){
						
						insert_node (best_ip, best_jp, LOOP);
					}
					else
					{
						fprintf (stderr, "NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
			//else if (type == MULTI)
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					int k, best_k = -1, best_row = -1;
					int tmp= INF, min = INF;
					for (k = i+1; k <= j-1; k++)
					{
						tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-1);
						if (tmp < min)
						{
							min = tmp;
							best_k = k;
							best_row = 1;
						}
						// TODO:
						// Hosna, May 1st, 2012
						// do I need to check for non-canonical base pairings here as well so the dangle values not be INF??
						if (fres[i+1].pair <= -1)
						{
							//AP
							for (auto &model : *energy_models) {
								model.energy_value = vm->get_energy_WM (i+2,k) +
									vm->get_energy_WM (k+1, j-1) +
									model.dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
									model.misc.multi_free_base_penalty;
							}
							tmp = emodel_energy_function (i, j, energy_models);

							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 2;
							}
						}
						if (fres[j-1].pair <= -1)
						{
							//AP
							for (auto &model : *energy_models) {
								model.energy_value = vm->get_energy_WM (i+1,k) +
									vm->get_energy_WM (k+1, j-2) +
									model.dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
									model.misc.multi_free_base_penalty;
							}
							tmp = emodel_energy_function (i, j, energy_models);

							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 3;
							}
						}
						if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
						{
							//AP
							for (auto &model : *energy_models) {
								model.energy_value = vm->get_energy_WM (i+2,k) +
									vm->get_energy_WM (k+1, j-2) +
									model.dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
									model.dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
									2*model.misc.multi_free_base_penalty;
							}
							tmp = emodel_energy_function (i, j, energy_models);

							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 4;
							}
						}

						// Hosna: June 28, 2007
						// the last branch of VM, which is WMB_(i+1),(j-1)
						// Hosna: July 5th, 2007
						// add a b_penalty to this case to match the previous cases
						tmp = WMB->get_energy(i+1,j-1)+ a_penalty + PSM_penalty+ b_penalty;
						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
						}

					}
					//printf("loop best row: %d, i:%d j:%d\n",best_row,i,j);
					switch (best_row)
					{
					case 1:
		        //      	printf("M_WM(%d,%d) branch 1: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-1);
						insert_node (i+1, best_k, M_WM);
						insert_node (best_k+1, j-1, M_WM);
						break;
					case 2:
		  //            	printf("M_WM(%d,%d) branch 2: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-1);
						insert_node (i+2, best_k, M_WM);
						insert_node (best_k+1, j-1, M_WM);
						break;
					case 3:
		   //           	printf("M_WM(%d,%d) branch 3: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-2);
						insert_node (i+1, best_k, M_WM);
						insert_node (best_k+1, j-2, M_WM);
						break;
					case 4:
		//              	printf("M_WM(%d,%d) branch 4: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-2);
						insert_node (i+2, best_k, M_WM);
						insert_node (best_k+1, j-2, M_WM);
						break;
					// Hosna: June 28, 2007
					// the last branch of VM, which is WMB_(i+1),(j-1)
					case 5:
		  //            	printf("M_WM(%d,%d) branch 3: pushing P_WMB(%d,%d)\n", i,j,i+1,j-1);
						insert_node(i+1,j-1, P_WMB);
						break;
					default:
						fprintf(stderr, "ERROR backtrack loop has no best row\n");
						exit(10);
					}
				}
					break;
			}
		}
			break;
		case FREE:
		{
			int j = cur_interval->j;
			
			if (j==0) return;

			int min = INF, tmp, best_row = -1, i, best_i, acc, energy_ij;

			if (debug)
				printf ("\t(0,%d) FREE\n", j);
			//if(KEVIN_DEBUG)
			//	for (i=0; i < nb_nucleotides; i++){printf("%c i=%d f_type=%c f_pair=%d\n",restricted[i],i,fres[i].type, fres[i].pair);} //kevin debug
			// this case is for j unpaired, so I have to check that.
			if (fres[j].pair <= -1)
			{
				tmp = W[j-1];
				
				if (tmp < min)
				{
					min = tmp;
					best_row = 0;
				}
			}
			for (i=0; i<=j-1; i++)    // no TURN
			{
				//AP
				//if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
				//	continue;

				if (int_sequence[i] == X || int_sequence[i+1] == X || int_sequence[j] == X || int_sequence[j-1] == X)
                    continue;

				// Don't need to make sure i and j don't have to pair with something else
				//  it's INF, done in fold_sequence_restricted
				acc = (i-1>0) ? W[i-1] : 0;
				energy_ij = v->get_energy(i,j);
				
				if (energy_ij < INF)
				{
					//AP
					for (auto &emodel : *energy_models) {
						model = &emodel;
						model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i], int_sequence[j], model) + acc;
					}
					tmp = emodel_energy_function (i, j, energy_models);

					if (tmp < min)
					{
						min = tmp;
						best_i = i;
						best_row = 1;
					}
				}

				if (fres[i].pair <= -1)
				{
					energy_ij = v->get_energy(i+1,j);
					if (energy_ij < INF)
					{
						//AP
						for (auto &emodel : *energy_models) {
							model = &emodel;
							model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j], model) + acc;
							model->energy_value += model->dangle_bot [int_sequence[j]]
								[int_sequence[i+1]]
								[int_sequence[i]];
						}
						tmp = emodel_energy_function (i, j, energy_models);

						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 2;
						}
					}
				}
				if (fres[j].pair <= -1)
				{
					energy_ij = v->get_energy(i,j-1);
					
					if (energy_ij < INF)
					{
						
						//AP
						for (auto &emodel : *energy_models) {
							model = &emodel;
							model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i], int_sequence[j-1], model) + acc;
							model->energy_value += model->dangle_top [int_sequence[j-1]]
								[int_sequence[i]]
								[int_sequence[j]];
							
						}
						tmp = emodel_energy_function (i, j, energy_models);
						
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 3;
						}
					}
				}
				if (fres[i].pair <= -1 && fres[j].pair <= -1)
				{
					energy_ij = v->get_energy(i+1,j-1);
					
					if (energy_ij < INF)
					{
						//AP
						for (auto &emodel : *energy_models) {
							model = &emodel;
							model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j-1], model) + acc;
							model->energy_value += model->dangle_bot [int_sequence[j-1]]
								[int_sequence[i+1]]
								[int_sequence[i]];
							model->energy_value += model->dangle_top [int_sequence[j-1]]
								[int_sequence[i+1]]
								[int_sequence[j]];
						}
						tmp = emodel_energy_function (i, j, energy_models);
						
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 4;
						}
					}
				}
			}
			// Hosna: June 28, 2007
			// the last branch of W, which is WMB_i,j
	//        energy_ij = WMB->get_energy(0,j);
	//        if (energy_ij < INF){
	//          	tmp = energy_ij + PS_penalty;
	//           	if (tmp < min){
	//           		min = tmp;
	//           		best_row = 5;
	//           	}
	//        }
		// Hosna June 30, 2007
		// The following would not take care of when
		// we have some unpaired bases before the start of the WMB
		for (i=0; i<=j-1; i++)
		{
            if (int_sequence[i] == X || int_sequence[j] == X)
                    continue;

			// Hosna: July 9, 2007
			// We only chop W to W + WMB when the bases before WMB are free
			if (i == 0 || (WMB->is_weakly_closed(0,i-1) && WMB->is_weakly_closed(i,j))){

				acc = (i-1>0) ? W[i-1]: 0;

				energy_ij = WMB->get_energy(i,j);

				if (energy_ij < INF)
				{
					tmp = energy_ij + PS_penalty + acc;

					if (tmp < min)
					{
						min = tmp;
						best_row = 5;
						best_i = i;
					}
				}

				// I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
				if (fres[i].pair <= -1 && i+1 < j)
				{
					energy_ij = WMB->get_energy(i+1,j);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 6;
							best_i = i;
						}
					}
				}

				// I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
				if (fres[j].pair <= -1 && i < j-1)
				{
					energy_ij = WMB->get_energy(i,j-1);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 7;
							best_i = i;
						}
					}
				}

				if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
				{
					energy_ij = WMB->get_energy(i+1,j-1);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 8;
							best_i = i;
						}
					}
				}
			}
		}
		
			switch (best_row)
			{
				case 0:
					if(KEVIN_DEBUG)
						printf("W(%d) case 0: inserting Free (0,%d)\n",j,j-1);
					insert_node (0, j-1, FREE); break;
				case 1:
					if(KEVIN_DEBUG)
						printf("W(%d) case 1: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, LOOP);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 2:
					if(KEVIN_DEBUG)
						printf("W(%d) case 2: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, LOOP);
					if (best_i >= 0)// Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 3:
					if(KEVIN_DEBUG)
						printf("W(%d) case 3: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, LOOP);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 4:
					if(KEVIN_DEBUG)
						printf("W(%d) case 4: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
		
					insert_node (best_i+1, j-1, LOOP);
				
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				// Hosna: June 28, 2007
				// the last branch of W, which is WMB_i,j
				case 5:
					if(KEVIN_DEBUG)
						printf("W(%d) case 5: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, P_WMB);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 6:
					if(KEVIN_DEBUG)
						printf("W(%d) case 6: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 7:
					if(KEVIN_DEBUG)
						printf("W(%d) case 7: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, P_WMB);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 8:
					if(KEVIN_DEBUG)
						printf("W(%d) case 8: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				default:
                	fprintf(stderr, "ERROR backtrack free has no best row\n");
                	exit(10);
			}
		}
			break;
		case M_WM:
//  else if(cur_interval->type == M_WM)
		{
			
			int i = cur_interval->i;
			int j = cur_interval->j;
			int tmp, min = INF;
			int best_k = -1, best_row = -1;

			if (debug)
				printf ("\t (%d,%d) M_WM\n", i,j);


			//AP
			//if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
			//	return;

			//AP
			for (auto &emodel : *energy_models) {
				model = &emodel;
				model->energy_value = v->get_energy(i,j) +
					AU_penalty_emodel (int_sequence[i], int_sequence[j], model) +
					model->misc.multi_helix_penalty;
			}
			tmp = emodel_energy_function (i, j, energy_models);

			if (tmp < min)
			{
				min = tmp;
				best_row = 1;
			}
			if (fres[i].pair <= -1)
			{
				//AP
				for (auto &emodel : *energy_models) {
					model = &emodel;
					model->energy_value = v->get_energy(i+1,j) +
						AU_penalty_emodel (int_sequence[i+1], int_sequence[j], model) +
						model->misc.multi_helix_penalty +
						model->misc.multi_free_base_penalty;
					//Aug 17 2017 kevin and Mahyar
                	//modiefied the formula such that we only add dangle_bot when i,j,i+1 are not X to avoid seg fault
					if(int_sequence[i] != X && int_sequence[j] != X && int_sequence[i+1] != X){
						model->energy_value += model->dangle_bot [int_sequence[j]] [int_sequence[i+1]] [int_sequence[i]];
					}
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 2;
				}
			}
		
			if (fres[j].pair <= -1)
			{
				//AP
				for (auto &emodel : *energy_models) {
					model = &emodel;
					model->energy_value = v->get_energy(i,j-1) +
						AU_penalty_emodel (int_sequence[i], int_sequence[j-1], model) +
						model->misc.multi_helix_penalty +
						model->misc.multi_free_base_penalty;
					//Aug 17 2017 kevin and Mahyar
                	//modiefied the formula such that we only add dangle_top when i,j,j-1 are not X to avoid seg fault
					if(int_sequence[i] != X && int_sequence[j] != X && int_sequence[j-1] != X){
						model->energy_value += model->dangle_top [int_sequence[j-1]] [int_sequence[i]] [int_sequence[j]];
					}
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 3;
				}
			}
			
			if (fres[i].pair <= -1 && fres[j].pair <= -1)
			{
				//AP
				for (auto &emodel : *energy_models) {
					model = &emodel;
					model->energy_value = v->get_energy(i+1,j-1) +
						AU_penalty_emodel (int_sequence[i+1], int_sequence[j-1], model) +
						model->misc.multi_helix_penalty +
						2*model->misc.multi_free_base_penalty;
					//Aug 17 2017 kevin and Mahyar
                	//modiefied the formula such that we only add dangle_bot,dangle_top when i,j,i+1, j-1 are not X to avoid seg fault
					if(int_sequence[i] != X && int_sequence[j] != X && int_sequence[i+1] != X && int_sequence[j-1] != X){
						model->energy_value += model->dangle_bot [int_sequence[j-1]] [int_sequence[i+1]] [int_sequence[i]] +
							model->dangle_top [int_sequence[j-1]] [int_sequence[i+1]] [int_sequence[j]];
					}
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 4;
				}
			}
			if (fres[i].pair <= -1)
			{
				//AP
				for (auto &model : *energy_models) {
					model.energy_value = vm->get_energy_WM (i+1,j) + model.misc.multi_free_base_penalty;
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 5;
				}
			}
			if (fres[j].pair <= -1)
			{
				//AP
				for (auto &model : *energy_models) {
					model.energy_value = vm->get_energy_WM (i,j-1) + model.misc.multi_free_base_penalty;
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 6;
				}
			}

			for (int k=i; k < j; k++)
				{
					tmp = vm->get_energy_WM (i, k) + vm->get_energy_WM (k+1, j);
					if (tmp < min)
					{
						min = tmp;
						best_k = k;
						best_row = 7;
					}
				}
			// Hosna: June 28, 2007
			// the last branch of WW, which is WMB_i,j
			tmp = WMB->get_energy(i,j)+PSM_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 8;
			}

			//18 Aug 2017 kevin and Mahyar
			//added this to check if j is a X
			//if it is, move j till it is the first X so the next iteration can handle it properly
			int new_j = j;
			if(int_sequence[new_j] == X){
				while(int_sequence[new_j-1] == X){
					new_j--;
					best_row = 9;
				}    
				
			}

			//18 Aug 2017 kevin and Mahyar
			//added this to check if i is a X
			//if it is, move i till it is the last X so the next iteration can handle it properly
			int new_i = i;
			if(int_sequence[new_i] == X){
				while(int_sequence[new_i+1] == X){
					new_i++;
					best_row = 9;
				}    
			}
			//18 Aug 2017 kevin and Mahyar
			//added this to make sure after we move the new i and j, we do not cross
			//if it is, error
			if(new_i >= new_j){
				best_row = -1; //error
			}
//printf("M_WM best row: %d, new_i: %d new_j: %d\n",best_row,new_i,new_j);
			switch (best_row)
				{
				case 1: insert_node (i, j, LOOP); break;
				case 2: insert_node (i+1, j, LOOP); break;
				case 3: insert_node (i, j-1, LOOP); break;
				case 4: insert_node (i+1, j-1, LOOP); break;
				case 5:
					if (j-i-1 > 0)
					insert_node (i+1, j, M_WM);
					break;
				case 6:
					if (j-1-i > 0)
					insert_node (i, j-1, M_WM);
					break;
				case 7:
					if (best_k-i > 0)
					insert_node (i, best_k, M_WM);
					if (j-best_k-1 > 0)
					insert_node (best_k+1, j, M_WM);
					break;
				// Hosna: June 28, 2007
				// the last branch of W, which is WMB_i,j
				case 8:
					insert_node(i,j,P_WMB);
					break;
				case 9: //18 Aug 2017 kevin and Mahyar added this case to jump i or j when encounter X
					insert_node(new_i,new_j,M_WM);
					break;
				default:
					printf("i= %d j=%d\n",i,j);
					fprintf(stderr, "ERROR backtrack M_WM has no best row\n");
					exit(10);
				}
			}
			break;
	// Hosna: Feb 19th 2007
		case P_WMB:
   // else if(cur_interval->type == P_WMB)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WMBP:
	// Hosna: April 18th, 2007
	// changed WMB to case 2 and WMBP
   // else if(cur_interval->type == P_WMBP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VP:
	//else if(cur_interval->type == P_VP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VPP:
	//else if(cur_interval->type == P_VPP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WI:
	//else if(cur_interval->type == P_WI)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_BE:
	//else if(cur_interval->type == P_BE)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WIP:
	//else if(cur_interval->type == P_WIP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		default:
			printf("Should not be here!\n"); //This should be removed or changed as it does nothing useful.
			exit(10);
	}

}

//AP
void W_final::backtrack_restricted_pkonly_emodel (seq_interval *cur_interval, str_features *fres){
	char type;
	energy_model *model;

	//Hosna, March 8, 2012
	// changing nested if to switch for optimality
	switch (cur_interval->type) {
			// TODO:
			// April 3, 2012
			// for the pk only case, I don't think I need to change any part of the LOOP case
			// April 18, 2012
			// I think the closing base pairs of the loops are set previously, but the pkonly condition needs to be checked
			// for the inner base pairs.
		case LOOP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;

			if (i >= j)
				return;
			f[i].pair = j;
			f[j].pair = i;

			//AP
			//if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
			//	return;

			// Hosna Jun. 28 2007
			// if the pairing is part of original structure, put '(' and ')' in the structure
			// otherwise make it '[' and ']'
			if (fres[i].pair == j){
				structure[i] = '(';
				structure[j] = ')';
			}else{
				fprintf(stderr, "Base pairing between %d and %d is pseudoknot free and should not happen here!! \n",i,j);
				exit(0);
				//structure[i] = '[';
				//structure[j] = ']';
			}

			type = v->get_type (i,j);
			if (debug)
				printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
			// Hosna, March 8, 2012
			// changing nested ifs to switch for optimality
			switch (type){
				case STACK:
					//if (type == STACK)
				{
					f[i].type = STACK;
					f[j].type = STACK;
					if (i+1 < j-1 && fres[i+1].pair == j-1) //check for pkonly
						this->insert_node(i+1,j-1, LOOP);
					//insert_node (i+1, j-1, LOOP);
					else
					{
						fprintf (stderr, "pkonly NOT GOOD RESTR STACK (i+1 and j-1 are not paired!), i=%d, j=%d\n", i, j);
						exit (0);
					}

				}
					break;
				case HAIRP:
					//else if (type == HAIRP)
				{
					f[i].type = HAIRP;
					f[j].type = HAIRP;
				}
					break;
				case INTER:
					//else if (type == INTER)
				{
					f[i].type = INTER;
					f[j].type = INTER;
					// detect the other closing pair
					int ip, jp, best_ip, best_jp, minq;
					int tmp, min = INF;
					// Hosna, August 31, 2012
					// The following restriction misses the long restricted loops, so I am chaning it
					//for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1) ; ip++)  // the -TURN shouldn't be there

/*
					//todo kevin confirm
					//23 Aug 2017 kevin and Mahyar
					//variable to store linker_length such that we change the range of ip and jp so we can check if i,ip and jp,j is larger than MAXLOOP (aka 30) properly for the cases where X is between i,j
					//if X is between i,j we treat it as if X does not exist 
					int skip = 0;
									
					if(is_cross_model(i,j)){
						skip = linker_length;
				
					}
					
					//todo kevin confirm
					//23 Aug 2017 kevin and Mahyar
					//added +skip to i+MAXLOOP+1
					for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1+skip); ip++)

*/
					for (ip = i+1; ip <= j-2 ; ip++)  // the -TURN shouldn't be there
					{
						// Hosna, August 28, 2012
						// TODO: cannot understand why we have th efollowing calculations, as it makes the following case be missed!
						// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
						// (____(_____((__(_______)__))_______________________________)____)
						// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
						// in this example int(5,59,11,27) is falsely missed and is equal to INF
						// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
						//minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
/*
						//todo kevin confirm
						//23 Aug 2017 kevin and Mayahr
						//changed minq to new one (copied from simfold backtrack emodel)
						minq = MAX (j-i+ip-MAXLOOP-2-skip, ip+1);
*/
						minq = ip+1;
						for (jp = minq; jp < j; jp++)
						{
							
							if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
								continue;
								

							//tmp = VBI->get_energy_str (i,j,ip,jp);
							// Hosna, March 26, 2012
							// modified to accommodate non-canonical base pairing in restricted structure
							// April 18, 2012
							// added a condition for pkonly
							//AP
							for (auto &model : *energy_models) {
								//model.energy_value = (fres[ip].pair == jp && fres[jp].pair == ip)? VBI->get_energy_str_restricted_emodel (i, j, ip, jp, fres, &model);:INF;
								//17 Aug 2017 kevin and Mahyar
								//changed the above function call to this one so we dont re-caculate it and just look up the value
								model.energy_value = (fres[ip].pair == jp && fres[jp].pair == ip)? v->get_energy(ip,jp):INF;
							}
							tmp = emodel_energy_function (i, j, energy_models);

							if (tmp < min)
							{
								min = tmp;
								best_ip = ip;
								best_jp = jp;
							}
						}
					}
					if (best_ip < best_jp)
						insert_node (best_ip, best_jp, LOOP);
					else
					{
						fprintf (stderr, "pkonly NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
					//else if (type == MULTI)
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					int k, best_k = -1, best_row = -1;
					int tmp= INF, min = INF;
					for (k = i+1; k <= j-1; k++)
					{
						tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-1);
						if (tmp < min)
						{
							min = tmp;
							best_k = k;
							best_row = 1;
						}
						if (fres[i+1].pair <= -1)
						{
							//AP
							for (auto &model : *energy_models) {
								model.energy_value = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-1) +
								model.dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
								model.misc.multi_free_base_penalty;
							}
							tmp = emodel_energy_function (i, j, energy_models);

							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 2;
							}
						}
						if (fres[j-1].pair <= -1)
						{
							//AP
							for (auto &model : *energy_models) {
								model.energy_value = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-2) +
								model.dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
								model.misc.multi_free_base_penalty;
							}
							tmp = emodel_energy_function (i, j, energy_models);

							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 3;
							}
						}
						if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
						{
							//AP
							for (auto &model : *energy_models) {
								model.energy_value = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-2) +
								model.dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
								model.dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
								2*model.misc.multi_free_base_penalty;
							}
							tmp = emodel_energy_function (i, j, energy_models);

							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 4;
							}
						}

						// Hosna: June 28, 2007
						// the last branch of VM, which is WMB_(i+1),(j-1)
						// Hosna: July 5th, 2007
						// add a b_penalty to this case to match the previous cases
						tmp = WMB->get_energy(i+1,j-1)+ a_penalty + PSM_penalty+ b_penalty;
						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
						}

					}
					switch (best_row)
					{
						case 1:
							//              	printf("M_WM(%d,%d) branch 1: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-1);
							insert_node (i+1, best_k, M_WM);
							insert_node (best_k+1, j-1, M_WM);
							break;
						case 2:
							//              	printf("M_WM(%d,%d) branch 2: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-1);
							insert_node (i+2, best_k, M_WM);
							insert_node (best_k+1, j-1, M_WM);
							break;
						case 3:
							//              	printf("M_WM(%d,%d) branch 3: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-2);
							insert_node (i+1, best_k, M_WM);
							insert_node (best_k+1, j-2, M_WM);
							break;
						case 4:
							//              	printf("M_WM(%d,%d) branch 4: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-2);
							insert_node (i+2, best_k, M_WM);
							insert_node (best_k+1, j-2, M_WM);
							break;
							// Hosna: June 28, 2007
							// the last branch of VM, which is WMB_(i+1),(j-1)
						case 5:
							//              	printf("M_WM(%d,%d) branch 3: pushing P_WMB(%d,%d)\n", i,j,i+1,j-1);
							insert_node(i+1,j-1, P_WMB);
							break;
						default:
							fprintf(stderr, "ERROR pkonly backtrack loop has no best row\n");
							exit(10);
					}
				}
					break;
			}
		}
			break;
		case FREE:
		{
			int j = cur_interval->j;

			if (j==0) return;

			int min = INF, tmp, best_row = -1, i, best_i, acc, energy_ij;

			if (debug)
				printf ("\t(0,%d) FREE\n", j);

			// this case is for j unpaired, so I have to check that.
			if (fres[j].pair <= -1)
			{
				tmp = W[j-1];
				if (tmp < min)
				{
					min = tmp;
					best_row = 0;
				}
			}
			for (i=0; i<=j-1; i++)    // no TURN
			{
				//AP
				//if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
				//	continue;

				if (int_sequence[i] == X || int_sequence[i+1] == X || int_sequence[j] == X || int_sequence[j-1] == X)
                    continue;

				// Don't need to make sure i and j don't have to pair with something else
				//  it's INF, done in fold_sequence_restricted
				acc = (i-1>0) ? W[i-1] : 0;
				// pkonly condition:
				energy_ij = (fres[j].pair == i && fres[i].pair== j)? v->get_energy(i,j): INF;
			
				if (energy_ij < INF)
				{
					//AP
					for (auto &emodel : *energy_models) {
						model = &emodel;
						model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i], int_sequence[j], model) + acc;
					}
					tmp = emodel_energy_function (i, j, energy_models);

					if (tmp < min)
					{
						min = tmp;
						best_i = i;
						best_row = 1;
					}
				}

				if (fres[i].pair <= -1)
				{
					energy_ij = (fres[j].pair == i+1 && fres[i+1].pair== j) ? v->get_energy(i+1,j) : INF;
					if (energy_ij < INF)
					{
						//AP
						for (auto &emodel : *energy_models) {
							model = &emodel;
							model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j], model) + acc;
							model->energy_value += model->dangle_bot [int_sequence[j]]
																	[int_sequence[i+1]]
																	[int_sequence[i]];
						}
						tmp = emodel_energy_function (i, j, energy_models);

						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 2;
						}
					}
				}
				if (fres[j].pair <= -1)
				{
					energy_ij = (fres[j-1].pair == i && fres[i].pair== j-1)? v->get_energy(i,j-1) : INF;
					if (energy_ij < INF)
					{
						//AP
						for (auto &emodel : *energy_models) {
							model = &emodel;
							model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i], int_sequence[j-1], model) + acc;
							model->energy_value += model->dangle_top [int_sequence[j-1]]
																	[int_sequence[i]]
																	[int_sequence[j]];
						}
						tmp = emodel_energy_function (i, j, energy_models);

						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 3;
						}
					}
				}
				if (fres[i].pair <= -1 && fres[j].pair <= -1)
				{
					energy_ij = (fres[j-1].pair == i+1 && fres[i+1].pair == j-1) ? v->get_energy(i+1,j-1) : INF;
					if (energy_ij < INF)
					{
						//AP
						for (auto &emodel : *energy_models) {
							model = &emodel;
							model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j-1], model) + acc;
							model->energy_value += model->dangle_bot [int_sequence[j-1]]
																	[int_sequence[i+1]]
																	[int_sequence[i]];
							model->energy_value += model->dangle_top [int_sequence[j-1]]
																	[int_sequence[i+1]]
																	[int_sequence[j]];
						}
						tmp = emodel_energy_function (i, j, energy_models);

						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 4;
						}
					}
				}
			}

			// Hosna: June 28, 2007
			// the last branch of W, which is WMB_i,j
			//        energy_ij = WMB->get_energy(0,j);
			//        if (energy_ij < INF){
			//          	tmp = energy_ij + PS_penalty;
			//           	if (tmp < min){
			//           		min = tmp;
			//           		best_row = 5;
			//           	}
			//        }
			// Hosna June 30, 2007
			// The following would not take care of when
			// we have some unpaired bases before the start of the WMB
			for (i=0; i<=j-1; i++)
			{
				// Hosna: July 9, 2007
				// We only chop W to W + WMB when the bases before WMB are free
				if (i == 0 || (WMB->is_weakly_closed(0,i-1) && WMB->is_weakly_closed(i,j))){

					acc = (i-1>0) ? W[i-1]: 0;

					energy_ij = WMB->get_energy(i,j);

					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;

						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
							best_i = i;
						}
					}

					// I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
					if (fres[i].pair <= -1 && i+1 < j)
					{
						energy_ij = WMB->get_energy(i+1,j);
						if (energy_ij < INF)
						{
							tmp = energy_ij + PS_penalty + acc;
							if (tmp < min)
							{
								min = tmp;
								best_row = 6;
								best_i = i;
							}
						}
					}

					// I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
					if (fres[j].pair <= -1 && i < j-1)
					{
						energy_ij = WMB->get_energy(i,j-1);
						if (energy_ij < INF)
						{
							tmp = energy_ij + PS_penalty + acc;
							if (tmp < min)
							{
								min = tmp;
								best_row = 7;
								best_i = i;
							}
						}
					}

					if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
					{
						energy_ij = WMB->get_energy(i+1,j-1);
						if (energy_ij < INF)
						{
							tmp = energy_ij + PS_penalty + acc;
							if (tmp < min)
							{
								min = tmp;
								best_row = 8;
								best_i = i;
							}
						}
					}
				}
			}
			switch (best_row)
			{
				case 0:
					//printf("W(%d) case 0: inserting Free (0,%d)\n",j,j-1);
					insert_node (0, j-1, FREE); break;
				case 1:
					//printf("W(%d) case 1: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, LOOP);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 2:
					//printf("W(%d) case 2: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, LOOP);
					if (best_i >= 0)// Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 3:
					//printf("W(%d) case 3: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, LOOP);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 4:
					//printf("W(%d) case 4: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, LOOP);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
					// Hosna: June 28, 2007
					// the last branch of W, which is WMB_i,j
				case 5:
					//printf("W(%d) case 5: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, P_WMB);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 6:
					//printf("W(%d) case 6: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 7:
					//printf("W(%d) case 7: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, P_WMB);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 8:
					//printf("W(%d) case 8: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				default:
					fprintf(stderr, "ERROR pkonly backtrack free has no best row\n");
					exit(10);
			}
		}
			break;
		case M_WM:
			//  else if(cur_interval->type == M_WM)
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			int tmp, min = INF;
			int best_k, best_row = -1;

			if (debug)
				printf ("\t (%d,%d) M_WM\n", i,j);

            // AP
			//if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
			//	return;

			// the if clause added April 3, 2012
			if (fres[j].pair == i && fres[i].pair == j){
				//AP
				for (auto &emodel : *energy_models) {
					model = &emodel;
					model->energy_value = v->get_energy(i,j) +
						AU_penalty_emodel (int_sequence[i], int_sequence[j], model) +
						model->misc.multi_helix_penalty;
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 1;
				}
			}
			if (fres[i].pair <= -1 && fres[j].pair == i+1 && fres[i+1].pair == j)
			{
				//AP
				for (auto &emodel : *energy_models) {
					model = &emodel;
					model->energy_value = v->get_energy(i+1,j) +
						AU_penalty_emodel (int_sequence[i+1], int_sequence[j], model) +
						model->misc.multi_helix_penalty +
						model->misc.multi_free_base_penalty;
					//Aug 18 2017 kevin and Mahyar
                	//modiefied the formula such that we only add dangle_bot when i,j,i+1 are not X to avoid seg fault
					if(int_sequence[i] != X && int_sequence[j] != X && int_sequence[i+1] != X){
						model->energy_value += model->dangle_bot [int_sequence[j]] [int_sequence[i+1]] [int_sequence[i]];
					}
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 2;
				}
			}
			if (fres[j].pair <= -1 && fres[j-1].pair == i && fres[i].pair == j-1)
			{
				//AP
				for (auto &emodel : *energy_models) {
					model = &emodel;
					model->energy_value = v->get_energy(i,j-1) +
						AU_penalty_emodel (int_sequence[i], int_sequence[j-1], model) +
						model->misc.multi_helix_penalty +
						model->misc.multi_free_base_penalty;
					//Aug 18 2017 kevin and Mahyar
                	//modiefied the formula such that we only add dangle_top when i,j,j-1 are not X to avoid seg fault
					if(int_sequence[i] != X && int_sequence[j] != X && int_sequence[j-1] != X){
						model->energy_value += model->dangle_top [int_sequence[j-1]] [int_sequence[i]] [int_sequence[j]];
					}
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 3;
				}
			}
			if (fres[i].pair <= -1 && fres[j].pair <= -1 && fres[j-1].pair == i+1 && fres[i+1].pair == j-1)
			{
				//AP
				for (auto &emodel : *energy_models) {
					model = &emodel;
					model->energy_value = v->get_energy(i+1,j-1) +
						AU_penalty_emodel (int_sequence[i+1], int_sequence[j-1], model) +
						model->misc.multi_helix_penalty +
						2*model->misc.multi_free_base_penalty;
					//Aug 18 2017 kevin and Mahyar
                	//modiefied the formula such that we only add dangle_bot,dangle_top when i,j,i+1, j-1 are not X to avoid seg fault
					if(int_sequence[i] != X && int_sequence[j] != X && int_sequence[i+1] != X && int_sequence[j-1] != X){
						model->energy_value += model->dangle_bot [int_sequence[j-1]] [int_sequence[i+1]] [int_sequence[i]] +
							model->dangle_top [int_sequence[j-1]] [int_sequence[i+1]] [int_sequence[j]];
					}
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 4;
				}
			}

			// TODO: April 3, 2012
			// do I need to change WM to pk_only as well?
			if (fres[i].pair <= -1)
			{
				//AP
				for (auto &model : *energy_models) {
					model.energy_value = vm->get_energy_WM (i+1,j) + model.misc.multi_free_base_penalty;
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 5;
				}
			}
			if (fres[j].pair <= -1)
			{
				//AP
				for (auto &model : *energy_models) {
					model.energy_value = vm->get_energy_WM (i,j-1) + model.misc.multi_free_base_penalty;
				}
				tmp = emodel_energy_function (i, j, energy_models);

				if (tmp < min)
				{
					min = tmp;
					best_row = 6;
				}
			}

			for (int k=i; k < j; k++)
			{
				tmp = vm->get_energy_WM (i, k) + vm->get_energy_WM (k+1, j);
				if (tmp < min)
				{
					min = tmp;
					best_k = k;
					best_row = 7;
				}
			}
			// Hosna: June 28, 2007
			// the last branch of W, which is WMB_i,j
			tmp = WMB->get_energy(i,j)+PSM_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 8;
			}

			//18 Aug 2017 kevin and Mahyar
			//added this to check if i is a X
			//if it is, move i till it is the last X so the next iteration can handle it properly
			int new_j = j;
			if(int_sequence[new_j] == X){
				while(int_sequence[new_j-1] == X){
					new_j--;
					best_row = 9;
				}    
				
			}

			//18 Aug 2017 kevin and Mahyar
			//added this to make sure after we move the new i and j, we do not cross
			//if it is, error
			int new_i = i;
			if(int_sequence[new_i] == X){
				while(int_sequence[new_i+1] == X){
					new_i++;
					best_row = 9;
				}    
			}

			//18 Aug 2017 kevin and Mahyar
			//added this to make sure after we move the new i and j, we do not cross
			//if it is, error
			if(new_i >= new_j){
				best_row = -1; //error
			}

			switch (best_row)
			{
				case 1: insert_node (i, j, LOOP); break;
				case 2: insert_node (i+1, j, LOOP); break;
				case 3: insert_node (i, j-1, LOOP); break;
				case 4: insert_node (i+1, j-1, LOOP); break;
				case 5:
					if (j-i-1 > 0)
						insert_node (i+1, j, M_WM);
					break;
				case 6:
					if (j-1-i > 0)
						insert_node (i, j-1, M_WM);
					break;
				case 7:
					if (best_k-i > 0)
						insert_node (i, best_k, M_WM);
					if (j-best_k-1 > 0)
						insert_node (best_k+1, j, M_WM);
					break;
					// Hosna: June 28, 2007
					// the last branch of W, which is WMB_i,j
				case 8:
					insert_node(i,j,P_WMB);
					break;
				case 9: //18 Aug 2017 kevin and Mahyar added this case to jump i or j when encounter X
					insert_node(new_i,new_j,M_WM);
					break;
				default:
					fprintf(stderr, "ERROR pkonly backtrack M_WM has no best row\n");
					exit(10);
			}
		}
			break;
			// Hosna: Feb 19th 2007
		case P_WMB:
			// else if(cur_interval->type == P_WMB)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WMBP:
			// Hosna: April 18th, 2007
			// changed WMB to case 2 and WMBP
			// else if(cur_interval->type == P_WMBP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VP:
			//else if(cur_interval->type == P_VP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VPP:
			//else if(cur_interval->type == P_VPP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WI:
			//else if(cur_interval->type == P_WI)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_BE:
			//else if(cur_interval->type == P_BE)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WIP:
			//else if(cur_interval->type == P_WIP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly_emodel(structure,f,cur_interval,energy_models);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		default:
			fprintf(stderr, "pkonly Should not be here!\n");
			exit(10);
	}
}

//Kevin and Ian July 24
void W_final::backtrack_restricted_simfold_emodel (seq_interval *cur_interval, str_features *fres)
// PRE:  All matrixes V, VM, WM and W have been filled
// POST: Discover the MFE path
{
    char type;
	energy_model *model;
    if(cur_interval->type == LOOP)
    {
        int i = cur_interval->i;
        int j = cur_interval->j;
        if (i >= j)
            return;
        f[i].pair = j;
        f[j].pair = i;
        structure[i] = '(';
        structure[j] = ')';

        type = V->get_type (i,j);
        if (debug)
            printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
        if (type == STACK)
        {
            f[i].type = STACK;
            f[j].type = STACK;
            if (i+1 < j-1)
                insert_node (i+1, j-1, LOOP);
            else
            {
                fprintf (stderr, "NOT GOOD RESTR STACK, i=%d, j=%d\n", i, j);
                exit (0);
            }

        }
        else if (type == HAIRP)
        {
            f[i].type = HAIRP;
            f[j].type = HAIRP;
        }
        else if (type == INTER)
        {
            f[i].type = INTER;
            f[j].type = INTER;
            // detect the other closing pair
            int ip, jp, best_ip = -1, best_jp = -1, minq;
            PARAMTYPE tmp, min = INF;

			//todo kevin confirm
			//23 Aug 2017 kevin and Mahyar
			//variable to store linker_length such that we change the range of ip and jp so we can check if i,ip and jp,j is larger than MAXLOOP (aka 30) properly for the cases where X is between i,j
			//if X is between i,j we treat it as if X does not exist 
			int skip = 0;
			
			if(is_cross_model(i,j)){
				skip = linker_length;
				
			}
			
			//todo kevin confirm
			//23 Aug 2017 kevin and Mahyar
			//added +skip to i+MAXLOOP+1
            for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1+skip); ip++)
            {
				//todo kevin confirm
				//23 Aug 2017 kevin and Mahyar
				//added -skip to j-i+ip-MAXLOOP-2
                minq = MAX (j-i+ip-MAXLOOP-2-skip, ip+1);

                for (jp = minq; jp < j; jp++)
                {

                    if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
                        continue;
					
                    //AP
                    for (auto &model : *energy_models) {
						//VBI->get_energy_str_restricted_emodel (i, j, ip, jp, fres, &model);
						//17 Aug 2017 kevin and Mahyar
						//changed the above function call to this one so we dont re-caculate it and just look up the value
                        model.energy_value = V->get_energy(ip,jp);  
						
                    }
                    tmp = emodel_energy_function (i, j, energy_models);

                    if (tmp < min)
                    {
                        min = tmp;
                        best_ip = ip;
                        best_jp = jp;
                    }
                }
            }
            if (best_ip < best_jp)
                insert_node (best_ip, best_jp, LOOP);
            else
            {
                fprintf (stderr, "NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
                exit (0);
            }
        }
        else if (type == MULTI)
        {
            f[i].type = MULTI;
            f[j].type = MULTI;
            int k, best_k = -1, best_row = -1;
            PARAMTYPE tmp, min = INF;
            for (k = i+1; k <= j-1; k++)
              {
                tmp = VM->get_energy_WM (i+1,k) + VM->get_energy_WM (k+1, j-1);
                if (tmp < min)
                  {
                    min = tmp;
                    best_k = k;
                    best_row = 1;
                  }
                if (fres[i+1].pair <= -1)
                {
                    //AP
                    for (auto &model : *energy_models) {
                        model.energy_value = VM->get_energy_WM (i+2,k) +
                            VM->get_energy_WM (k+1, j-1) +
                            model.dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
                            model.misc.multi_free_base_penalty;
                    }
                    tmp = emodel_energy_function (i, j, energy_models);

                    if (tmp < min)
                    {
                        min = tmp;
                        best_k = k;
                        best_row = 2;
                    }
                }
                if (fres[j-1].pair <= -1)
                {
                    //AP
                    for (auto &model : *energy_models) {
						model.energy_value = VM->get_energy_WM (i+1,k) +
							VM->get_energy_WM (k+1, j-2) +
							model.dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
							model.misc.multi_free_base_penalty;
					}
                    tmp = emodel_energy_function (i, j, energy_models);

                    if (tmp < min)
                    {
                        min = tmp;
                        best_k = k;
                        best_row = 3;
                    }
                }
                if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
                {
                    //AP
					for (auto &model : *energy_models) {
						model.energy_value = VM->get_energy_WM (i+2,k) + VM->get_energy_WM (k+1, j-2) +
                    		model.dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
                    		model.dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
                    		2*model.misc.multi_free_base_penalty;
                    }
							tmp = emodel_energy_function (i, j, energy_models);
					if (tmp < min)
                    {
                        min = tmp;
                        best_k = k;
                        best_row = 4;
                    }
                }
              }
            switch (best_row)
              {
              case 1: insert_node (i+1, best_k, M_WM);
                insert_node (best_k+1, j-1, M_WM); break;
              case 2: insert_node (i+2, best_k, M_WM);
                insert_node (best_k+1, j-1, M_WM); break;
              case 3: insert_node (i+1, best_k, M_WM);
                insert_node (best_k+1, j-2, M_WM); break;
              case 4: insert_node (i+2, best_k, M_WM);
                insert_node (best_k+1, j-2, M_WM); break;
			  default:
				fprintf(stderr, "ERROR backtrack loop has no best row\n");
				exit(10);
              }
        }
    }
    else if(cur_interval->type == FREE)
    {
        int j = cur_interval->j;

        if (j==0) return;
        //if (j <= TURN) return;

        PARAMTYPE min = INF, tmp, acc, energy_ij;
        int best_row = -1, best_i = -1;

        if (debug)
            printf ("\t(0,%d) FREE\n", j);

        // this case is for j unpaired, so I have to check that.
        if (fres[j].pair <= -1)
        {
        //printf ("j=%d\n", j);
            tmp = W[j-1];

	   // Ian Wark July 26 2017
           // Removed if statement.
           // If this is unpaired, tmp will be INF
           // if tmp is INF, it cannot be less than min, which is default INF
           // So if unpaired best_row will not be properly set to 0
           // and it will not recognized that it should be unpaired.
           // TODO this error will also be in normal simfold backtrack and should be changed
           //if (tmp < min)
           //{
                min = tmp;
                best_row = 0;
            //}
        }
        for (int i=0; i<=j-1; i++)    // no TURN
        {
			//AP
			//if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
			//	continue;

			if (int_sequence[i] == X || int_sequence[j] == X)
                continue;

            // Don't need to make sure i and j don't have to pair with something else
            //  it's INF, done in fold_sequence_restricted
            acc = (i-1>0) ? W[i-1] : 0;
            energy_ij = V->get_energy(i,j);
            if (energy_ij < INF)
            {
                //AP
				for (auto &emodel : *energy_models) {
					model = &emodel;
					model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i], int_sequence[j], model) + acc;
				}
				tmp = emodel_energy_function (i, j, energy_models);
//                 if (fres[i].pair == j)
//                 {
//                     min = tmp;
//                     best_i = i;
//                     best_row = 1;
//                     continue;
//                 }
                if (tmp < min)
                {
                min = tmp;
                best_i = i;
                best_row = 1;
                }
            }

            if (fres[i].pair <= -1)
            {
                energy_ij = V->get_energy(i+1,j);
                if (energy_ij < INF)
                {
                    //AP
					for (auto &emodel : *energy_models) {
						model = &emodel;
						model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j], model) + acc;
						model->energy_value += model->dangle_bot [int_sequence[j]]
							[int_sequence[i+1]]
							[int_sequence[i]];
					}
					tmp = emodel_energy_function (i, j, energy_models);
    //                 if (fres[i+1].pair == j)
    //                 {
    //                     min = tmp;
    //                     best_i = i;
    //                     best_row = 2;
    //                     continue;
    //                 }
                    if (tmp < min)
                    {
                        min = tmp;
                        best_i = i;
                        best_row = 2;
                    }
                }
            }
            if (fres[j].pair <= -1)
            {
                energy_ij = V->get_energy(i,j-1);
                if (energy_ij < INF)
                {
                    //AP
					for (auto &emodel : *energy_models) {
						model = &emodel;
						model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i], int_sequence[j-1], model) + acc;
						model->energy_value += model->dangle_top [int_sequence[j-1]]
							[int_sequence[i]]
							[int_sequence[j]];
					}
					tmp = emodel_energy_function (i, j, energy_models);
    //                 if (fres[i].pair == j-1)
    //                 {
    //                     min = tmp;
    //                     best_i = i;
    //                     best_row = 3;
    //                     continue;
    //                 }
                    if (tmp < min)
                    {
                        min = tmp;
                        best_i = i;
                        best_row = 3;
                    }
                }
            }
            if (fres[i].pair <= -1 && fres[j].pair <= -1)
            {
                energy_ij = V->get_energy(i+1,j-1);

                if (energy_ij < INF)
                {
                    //AP
					for (auto &emodel : *energy_models) {
						model = &emodel;
						model->energy_value = energy_ij + AU_penalty_emodel (int_sequence[i+1], int_sequence[j-1], model) + acc;
						model->energy_value += model->dangle_bot [int_sequence[j-1]]
							[int_sequence[i+1]]
							[int_sequence[i]];
						model->energy_value += model->dangle_top [int_sequence[j-1]]
							[int_sequence[i+1]]
							[int_sequence[j]];
					}
					tmp = emodel_energy_function (i, j, energy_models);
    //                 if (fres[i+1].pair == j-1)
    //                 {
    //                     min = tmp;
    //                     best_i = i;
    //                     best_row = 4;
    //                     continue;
    //                 }
                    if (tmp < min)
                    {
                        min = tmp;
                        best_i = i;
                        best_row = 4;
                    }
                }
            }
        }

        switch (best_row)
        {
            case 0: insert_node (0, j-1, FREE); break;
            case 1: insert_node (best_i, j, LOOP);
                if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
                    insert_node (0, best_i-1, FREE);
                break;
            case 2: insert_node (best_i+1, j, LOOP);
                if (best_i-1 > 0)
                    insert_node (0, best_i-1, FREE);
                break;
            case 3: insert_node (best_i, j-1, LOOP);
                if (best_i-1 > 0)
                    insert_node (0, best_i-1, FREE);
                break;
            case 4: insert_node (best_i+1, j-1, LOOP);
                if (best_i-1 > 0)
                    insert_node (0, best_i-1, FREE);
                break;
            default:
                fprintf(stderr, "ERROR backtrack free has no best row\n");
                exit(10);
        }
    }
  else if(cur_interval->type == M_WM)
    {
      int i = cur_interval->i;
      int j = cur_interval->j;
      PARAMTYPE tmp, min = INF;
      int best_k = -1, best_row = -1;

      if (debug)
        printf ("\t (%d,%d) M_WM\n", i,j);

        //AP
        //if (int_sequence[i] == X || int_sequence[j] == X || int_sequence[i+1] == X || int_sequence[j+1] == X || int_sequence[i-1] == X || int_sequence[j-1] == X)
        //	return;

	//AP
	for (auto &emodel : *energy_models) {
		model = &emodel;
		model->energy_value = V->get_energy(i,j) +
			AU_penalty_emodel (int_sequence[i], int_sequence[j], model) +
			model->misc.multi_helix_penalty;
	}
	tmp = emodel_energy_function (i, j, energy_models);
      if (tmp < min)
        {
          min = tmp;
          best_row = 1;
        }
      if (fres[i].pair <= -1)
      {
        //AP
		for (auto &emodel : *energy_models) {
			model = &emodel;
			model->energy_value = V->get_energy(i+1,j) +
				AU_penalty_emodel (int_sequence[i+1], int_sequence[j], model) +
				model->misc.multi_helix_penalty +
				model->misc.multi_free_base_penalty;
			//Aug 18 2017 kevin and Mahyar
			//modiefied the formula such that we only add dangle_bot when i,j,i+1 are not X to avoid seg fault
			if(int_sequence[i] != X && int_sequence[j] != X && int_sequence[i+1] != X){
				model->energy_value += model->dangle_bot [int_sequence[j]] [int_sequence[i+1]] [int_sequence[i]];
			}
		}
		tmp = emodel_energy_function (i, j, energy_models);
          if (tmp < min)
          {
              min = tmp;
              best_row = 2;
          }
      }
      if (fres[j].pair <= -1)
      {
        //AP
		for (auto &emodel : *energy_models) {
			model = &emodel;
			model->energy_value = V->get_energy(i,j-1) +
				AU_penalty_emodel (int_sequence[i], int_sequence[j-1], model) +
				model->misc.multi_helix_penalty +
				model->misc.multi_free_base_penalty;
			//Aug 18 2017 kevin and Mahyar
			//modiefied the formula such that we only add dangle_top when i,j,j-1 are not X to avoid seg fault
			if(int_sequence[i] != X && int_sequence[j] != X && int_sequence[j-1] != X){
				model->energy_value += model->dangle_top [int_sequence[j-1]] [int_sequence[i]] [int_sequence[j]];
			}
		}
		tmp = emodel_energy_function (i, j, energy_models);
          if (tmp < min)
          {
              min = tmp;
              best_row = 3;
          }
      }
      if (fres[i].pair <= -1 && fres[j].pair <= -1)
      {
          //AP
				for (auto &emodel : *energy_models) {
					model = &emodel;
					model->energy_value = V->get_energy(i+1,j-1) +
						AU_penalty_emodel (int_sequence[i+1], int_sequence[j-1], model) +
						model->misc.multi_helix_penalty +
						2*model->misc.multi_free_base_penalty;
					//Aug 18 2017 kevin and Mahyar
                	//modiefied the formula such that we only add dangle_bot,dangle_top when i,j,i+1, j-1 are not X to avoid seg fault
					if(int_sequence[i] != X && int_sequence[j] != X && int_sequence[i+1] != X && int_sequence[j-1] != X){
						model->energy_value += model->dangle_bot [int_sequence[j-1]] [int_sequence[i+1]] [int_sequence[i]] +
							model->dangle_top [int_sequence[j-1]] [int_sequence[i+1]] [int_sequence[j]];
					}
				}
				tmp = emodel_energy_function (i, j, energy_models);
          if (tmp < min)
          {
              min = tmp;
              best_row = 4;
          }
      }
      if (fres[i].pair <= -1)
      {
          //AP
		for (auto &model : *energy_models) {
			model.energy_value = VM->get_energy_WM (i+1,j) + model.misc.multi_free_base_penalty;
		}
		tmp = emodel_energy_function (i, j, energy_models);
          if (tmp < min)
          {
              min = tmp;
              best_row = 5;
          }
      }
      if (fres[j].pair <= -1)
      {
          //AP
		for (auto &model : *energy_models) {
			model.energy_value = VM->get_energy_WM (i,j-1) + model.misc.multi_free_base_penalty;
		}
		tmp = emodel_energy_function (i, j, energy_models);
          if (tmp < min)
          {
              min = tmp;
              best_row = 6;
          }
      }

      for (int k=i; k < j; k++)
        {
            tmp = VM->get_energy_WM (i, k) + VM->get_energy_WM (k+1, j);
            if (tmp < min)
              {
                min = tmp;
                best_k = k;
                best_row = 7;
              }
        }
		
		//18 Aug 2017 kevin and Mahyar
		//added this to check if j is a X
		//if it is, move j till it is the first X so the next iteration can handle it properly
		int new_j = j;
		if(int_sequence[new_j] == X){
			while(int_sequence[new_j-1] == X){
				new_j--;
				best_row = 9;
			}    
			
		}

		//18 Aug 2017 kevin and Mahyar
		//added this to check if i is a X
		//if it is, move i till it is the last X so the next iteration can handle it properly
		int new_i = i;
		if(int_sequence[new_i] == X){
			while(int_sequence[new_i+1] == X){
				new_i++;
				best_row = 9;
			}    
		}
		//18 Aug 2017 kevin and Mahyar
		//added this to make sure after we move the new i and j, we do not cross
		//if it is, error
		if(new_i >= new_j){
			best_row = -1; //error
		}

      switch (best_row)
        {
          case 1: insert_node (i, j, LOOP); break;
          case 2: insert_node (i+1, j, LOOP); break;
          case 3: insert_node (i, j-1, LOOP); break;
          case 4: insert_node (i+1, j-1, LOOP); break;
          case 5:
            if (j-i-1 > 0)
              insert_node (i+1, j, M_WM);
            break;
          case 6:
            if (j-1-i > 0)
              insert_node (i, j-1, M_WM);
            break;
          case 7:
            if (best_k-i > 0)
              insert_node (i, best_k, M_WM);
            if (j-best_k-1 > 0)
              insert_node (best_k+1, j, M_WM);
            break;
		  case 9: //18 Aug 2017 kevin and Mahyar added this case to jump i or j when encounter X
				insert_node(new_i,new_j,M_WM);
				break;
		  default:
			fprintf(stderr, "ERROR backtrack loop has no best row\n");
			exit(10);
          }
    }
}


void W_final::backtrack_restricted_pkonly (seq_interval *cur_interval, str_features *fres){
    char type;
	//Hosna, March 8, 2012
	// changing nested if to switch for optimality
	switch (cur_interval->type){
			// TODO:
			// April 3, 2012
			// for the pk only case, I don't think I need to change any part of the LOOP case
			// April 18, 2012
			// I think the closing base pairs of the loops are set previously, but the pkonly condition needs to be checked
			// for the inner base pairs.
		case LOOP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i >= j)
				return;
			f[i].pair = j;
			f[j].pair = i;

			// Hosna Jun. 28 2007
			// if the pairing is part of original structure, put '(' and ')' in the structure
			// otherwise make it '[' and ']'
			if (fres[i].pair == j){
				structure[i] = '(';
				structure[j] = ')';
			}else{
				fprintf(stderr, "Base pairing between %d and %d is pseudoknot free and should not happen here!! \n",i,j);
				exit(0);
				//structure[i] = '[';
				//structure[j] = ']';
			}

			type = v->get_type (i,j);
			if (debug)
				printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
			// Hosna, March 8, 2012
			// changing nested ifs to switch for optimality
			switch (type){
				case STACK:
					//if (type == STACK)
				{
					f[i].type = STACK;
					f[j].type = STACK;
					if (i+1 < j-1 && fres[i+1].pair == j-1) //check for pkonly
						this->insert_node(i+1,j-1, LOOP);
					//insert_node (i+1, j-1, LOOP);
					else
					{
						fprintf (stderr, "NOT GOOD RESTR STACK (i+1 and j-1 are not paired!), i=%d, j=%d\n", i, j);
						exit (0);
					}

				}
					break;
				case HAIRP:
					//else if (type == HAIRP)
				{
					f[i].type = HAIRP;
					f[j].type = HAIRP;
				}
					break;
				case INTER:
					//else if (type == INTER)
				{
					f[i].type = INTER;
					f[j].type = INTER;
					// detect the other closing pair
					int ip, jp, best_ip, best_jp, minq;
					int tmp, min = INF;
					// Hosna, August 31, 2012
					// The following restriction misses the long restricted loops, so I am chaning it
					//for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1) ; ip++)  // the -TURN shouldn't be there
					for (ip = i+1; ip <= j-2 ; ip++)  // the -TURN shouldn't be there
					{
						// Hosna, August 28, 2012
						// TODO: cannot understand why we have th efollowing calculations, as it makes the following case be missed!
						// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
						// (____(_____((__(_______)__))_______________________________)____)
						// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
						// in this example int(5,59,11,27) is falsely missed and is equal to INF
						// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
						//minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
						minq = ip+1;
						for (jp = minq; jp < j; jp++)
						{
							if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
								continue;
							//tmp = VBI->get_energy_str (i,j,ip,jp);
							// Hosna, March 26, 2012
							// modified to accommodate non-canonical base pairing in restricted structure
							// April 18, 2012
							// added a condition for pkonly
							tmp = (fres[ip].pair == jp && fres[jp].pair == ip)? VBI->get_energy_str_restricted(i,j,ip,jp,fres):INF;
							if (tmp < min)
							{
								min = tmp;
								best_ip = ip;
								best_jp = jp;
							}
						}
					}
					if (best_ip < best_jp)
						insert_node (best_ip, best_jp, LOOP);
					else
					{
						fprintf (stderr, "NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
					//else if (type == MULTI)
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					int k, best_k = -1, best_row = -1;
					int tmp= INF, min = INF;
					for (k = i+1; k <= j-1; k++)
					{
						tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-1);
						if (tmp < min)
						{
							min = tmp;
							best_k = k;
							best_row = 1;
						}
						if (fres[i+1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-1) +
							dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
							misc.multi_free_base_penalty;
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 2;
							}
						}
						if (fres[j-1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-2) +
							dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
							misc.multi_free_base_penalty;
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 3;
							}
						}
						if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-2) +
							dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
							dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
							2*misc.multi_free_base_penalty;
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 4;
							}
						}

						// Hosna: June 28, 2007
						// the last branch of VM, which is WMB_(i+1),(j-1)
						// Hosna: July 5th, 2007
						// add a b_penalty to this case to match the previous cases
						tmp = WMB->get_energy(i+1,j-1)+ a_penalty + PSM_penalty+ b_penalty;
						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
						}

					}
					switch (best_row)
					{
						case 1:
							//              	printf("M_WM(%d,%d) branch 1: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-1);
							insert_node (i+1, best_k, M_WM);
							insert_node (best_k+1, j-1, M_WM);
							break;
						case 2:
							//              	printf("M_WM(%d,%d) branch 2: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-1);
							insert_node (i+2, best_k, M_WM);
							insert_node (best_k+1, j-1, M_WM);
							break;
						case 3:
							//              	printf("M_WM(%d,%d) branch 3: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-2);
							insert_node (i+1, best_k, M_WM);
							insert_node (best_k+1, j-2, M_WM);
							break;
						case 4:
							//              	printf("M_WM(%d,%d) branch 4: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-2);
							insert_node (i+2, best_k, M_WM);
							insert_node (best_k+1, j-2, M_WM);
							break;
							// Hosna: June 28, 2007
							// the last branch of VM, which is WMB_(i+1),(j-1)
						case 5:
							//              	printf("M_WM(%d,%d) branch 3: pushing P_WMB(%d,%d)\n", i,j,i+1,j-1);
							insert_node(i+1,j-1, P_WMB);
							break;
					}
				}
					break;
			}
		}
			break;
		case FREE:
		{
			int j = cur_interval->j;

			if (j==0) return;

			int min = INF, tmp, best_row, i, best_i, acc, energy_ij;

			if (debug)
				printf ("\t(0,%d) FREE\n", j);

			// this case is for j unpaired, so I have to check that.
			if (fres[j].pair <= -1)
			{
				tmp = W[j-1];
				if (tmp < min)
				{
					min = tmp;
					best_row = 0;
				}
			}
			for (i=0; i<=j-1; i++)    // no TURN
			{

				// Don't need to make sure i and j don't have to pair with something else
				//  it's INF, done in fold_sequence_restricted
				acc = (i-1>0) ? W[i-1] : 0;
				// pkonly condition:
				energy_ij = (fres[j].pair == i && fres[i].pair== j)? v->get_energy(i,j): INF;
			
				if (energy_ij < INF)
				{
					tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
					if (tmp < min)
					{
						min = tmp;
						best_i = i;
						best_row = 1;
					}
				}

				if (fres[i].pair <= -1)
				{
					energy_ij = (fres[j].pair == i+1 && fres[i+1].pair== j) ? v->get_energy(i+1,j) : INF;
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
						tmp += dangle_bot [int_sequence[j]]
						[int_sequence[i+1]]
						[int_sequence[i]];
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 2;
						}
					}
				}
				if (fres[j].pair <= -1)
				{
					energy_ij = (fres[j-1].pair == i && fres[i].pair== j-1)? v->get_energy(i,j-1) : INF;
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j-1]) + acc;
						tmp += dangle_top [int_sequence[j-1]]
						[int_sequence[i]]
						[int_sequence[j]];
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 3;
						}
					}
				}
				if (fres[i].pair <= -1 && fres[j].pair <= -1)
				{
					energy_ij = (fres[j-1].pair == i+1 && fres[i+1].pair == j-1) ? v->get_energy(i+1,j-1) : INF;
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
						tmp += dangle_bot [int_sequence[j-1]]
						[int_sequence[i+1]]
						[int_sequence[i]];
						tmp += dangle_top [int_sequence[j-1]]
						[int_sequence[i+1]]
						[int_sequence[j]];
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 4;
						}
					}
				}
			}

			// Hosna: June 28, 2007
			// the last branch of W, which is WMB_i,j
			//        energy_ij = WMB->get_energy(0,j);
			//        if (energy_ij < INF){
			//          	tmp = energy_ij + PS_penalty;
			//           	if (tmp < min){
			//           		min = tmp;
			//           		best_row = 5;
			//           	}
			//        }
			// Hosna June 30, 2007
			// The following would not take care of when
			// we have some unpaired bases before the start of the WMB
			for (i=0; i<=j-1; i++)
			{
				// Hosna: July 9, 2007
				// We only chop W to W + WMB when the bases before WMB are free
				if (i == 0 || (WMB->is_weakly_closed(0,i-1) && WMB->is_weakly_closed(i,j))){

					acc = (i-1>0) ? W[i-1]: 0;

					energy_ij = WMB->get_energy(i,j);

					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;

						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
							best_i = i;
						}
					}

					// I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
					if (fres[i].pair <= -1 && i+1 < j)
					{
						energy_ij = WMB->get_energy(i+1,j);
						if (energy_ij < INF)
						{
							tmp = energy_ij + PS_penalty + acc;
							if (tmp < min)
							{
								min = tmp;
								best_row = 6;
								best_i = i;
							}
						}
					}

					// I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
					if (fres[j].pair <= -1 && i < j-1)
					{
						energy_ij = WMB->get_energy(i,j-1);
						if (energy_ij < INF)
						{
							tmp = energy_ij + PS_penalty + acc;
							if (tmp < min)
							{
								min = tmp;
								best_row = 7;
								best_i = i;
							}
						}
					}

					if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
					{
						energy_ij = WMB->get_energy(i+1,j-1);
						if (energy_ij < INF)
						{
							tmp = energy_ij + PS_penalty + acc;
							if (tmp < min)
							{
								min = tmp;
								best_row = 8;
								best_i = i;
							}
						}
					}
				}
			}
			switch (best_row)
			{
				case 0:
					//printf("W(%d) case 0: inserting Free (0,%d)\n",j,j-1);
					insert_node (0, j-1, FREE); break;
				case 1:
					//printf("W(%d) case 1: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, LOOP);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 2:
					//printf("W(%d) case 2: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, LOOP);
					if (best_i >= 0)// Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 3:
					//printf("W(%d) case 3: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, LOOP);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 4:
					//printf("W(%d) case 4: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, LOOP);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
					// Hosna: June 28, 2007
					// the last branch of W, which is WMB_i,j
				case 5:
					//printf("W(%d) case 5: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, P_WMB);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 6:
					//printf("W(%d) case 6: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 7:
					//printf("W(%d) case 7: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, P_WMB);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 8:
					//printf("W(%d) case 8: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
			}
		}
			break;
		case M_WM:
			//  else if(cur_interval->type == M_WM)
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			int tmp, min = INF;
			int best_k, best_row;

			if (debug)
				printf ("\t (%d,%d) M_WM\n", i,j);

			// the if clause added April 3, 2012
			if (fres[j].pair == i && fres[i].pair == j){
				tmp = v->get_energy(i,j) +
				AU_penalty (int_sequence[i], int_sequence[j]) +
				misc.multi_helix_penalty;
				if (tmp < min)
				{
					min = tmp;
					best_row = 1;
				}
			}
			if (fres[i].pair <= -1 && fres[j].pair == i+1 && fres[i+1].pair == j)
			{
				tmp = v->get_energy(i+1,j) +
				AU_penalty (int_sequence[i+1], int_sequence[j]) +
				dangle_bot [int_sequence[j]]
				[int_sequence[i+1]]
				[int_sequence[i]] +
				misc.multi_helix_penalty +
				misc.multi_free_base_penalty;
				if (tmp < min)
				{
					min = tmp;
					best_row = 2;
				}
			}
			if (fres[j].pair <= -1 && fres[j-1].pair == i && fres[i].pair == j-1)
			{
				tmp = v->get_energy(i,j-1) +
				AU_penalty (int_sequence[i], int_sequence[j-1]) +
				dangle_top [int_sequence[j-1]]
				[int_sequence[i]]
				[int_sequence[j]] +
				misc.multi_helix_penalty +
				misc.multi_free_base_penalty;
				if (tmp < min)
				{
					min = tmp;
					best_row = 3;
				}
			}
			if (fres[i].pair <= -1 && fres[j].pair <= -1 && fres[j-1].pair == i+1 && fres[i+1].pair == j-1)
			{
				tmp = v->get_energy(i+1,j-1) +
				AU_penalty (int_sequence[i+1], int_sequence[j-1]) +
				dangle_bot [int_sequence[j-1]]
				[int_sequence[i+1]]
				[int_sequence[i]] +
				dangle_top [int_sequence[j-1]]
				[int_sequence[i+1]]
				[int_sequence[j]] +
				misc.multi_helix_penalty +
				2*misc.multi_free_base_penalty;
				if (tmp < min)
				{
					min = tmp;
					best_row = 4;
				}
			}

			// TODO: April 3, 2012
			// do I need to change WM to pk_only as well?
			if (fres[i].pair <= -1)
			{
				tmp = vm->get_energy_WM (i+1,j) + misc.multi_free_base_penalty;
				if (tmp < min)
				{
					min = tmp;
					best_row = 5;
				}
			}
			if (fres[j].pair <= -1)
			{
				tmp = vm->get_energy_WM (i,j-1) + misc.multi_free_base_penalty;
				if (tmp < min)
				{
					min = tmp;
					best_row = 6;
				}
			}

			for (int k=i; k < j; k++)
			{
				tmp = vm->get_energy_WM (i, k) + vm->get_energy_WM (k+1, j);
				if (tmp < min)
				{
					min = tmp;
					best_k = k;
					best_row = 7;
				}
			}
			// Hosna: June 28, 2007
			// the last branch of W, which is WMB_i,j
			tmp = WMB->get_energy(i,j)+PSM_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 8;
			}

			switch (best_row)
			{
				case 1: insert_node (i, j, LOOP); break;
				case 2: insert_node (i+1, j, LOOP); break;
				case 3: insert_node (i, j-1, LOOP); break;
				case 4: insert_node (i+1, j-1, LOOP); break;
				case 5:
					if (j-i-1 > 0)
						insert_node (i+1, j, M_WM);
					break;
				case 6:
					if (j-1-i > 0)
						insert_node (i, j-1, M_WM);
					break;
				case 7:
					if (best_k-i > 0)
						insert_node (i, best_k, M_WM);
					if (j-best_k-1 > 0)
						insert_node (best_k+1, j, M_WM);
					break;
					// Hosna: June 28, 2007
					// the last branch of W, which is WMB_i,j
				case 8:
					insert_node(i,j,P_WMB);
					break;
			}
		}
			break;
			// Hosna: Feb 19th 2007
		case P_WMB:
			// else if(cur_interval->type == P_WMB)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WMBP:
			// Hosna: April 18th, 2007
			// changed WMB to case 2 and WMBP
			// else if(cur_interval->type == P_WMBP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VP:
			//else if(cur_interval->type == P_VP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VPP:
			//else if(cur_interval->type == P_VPP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WI:
			//else if(cur_interval->type == P_WI)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_BE:
			//else if(cur_interval->type == P_BE)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WIP:
			//else if(cur_interval->type == P_WIP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		default:
			printf("Should not be here!\n");
	}


}

void W_final::print_result ()
// PRE:  The matrix V has been calculated and the results written in f
// POST: Prints details of each elementary structure
{
    int i;
    int energy = INF, sum;

    printf ("Minimum energy: %d\n", W[nb_nucleotides-1]);
    sum = 0;

    for (i=0; i< nb_nucleotides; i++)
    {
        if (f[i].pair > i)
        {
			//Hosna March 8, 2012
			// changing nested ifs to switch for optimality
			switch (f[i].type){
				case HAIRP:
				//if (f[i].type == HAIRP)
					energy = V->get_energy(i, f[i].pair);
					break;
				case STACK:
				//else if (f[i].type == STACK)
					energy = V->get_energy(i, f[i].pair) - V->get_energy(i+1, f[i+1].pair);
					break;
				case P_VP:
				// Hosna: June 28th, 2007
				//else if (f[i].type == P_VP){
					energy = WMB->get_VP(i,f[i].pair);
					break;
				case P_VPP:
				//}else if(f[i].type == P_VPP){
					energy = WMB->get_VPP(i,f[i].pair);
				//}
					break;
			}
            printf ("Pair (%d,%d), type %c,\tenergy %6d\n", i, f[i].pair, f[i].type, energy);
            sum += energy;
        }
    }
    printf ("0....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8\n");
    printf ("%s\n", sequence);
    printf ("%s\n", structure);

}
