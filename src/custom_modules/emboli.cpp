/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./emboli.h"

// declare cell definitions here 

Cell_Definition epithelial; 
Cell_Definition mesenchymal;
Cell_Definition stem_like; 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 

	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it. 
	
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 

	int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index(PhysiCell_constants::G0G1_phase );
	int S_index = flow_cytometry_separated_cycle_model.find_phase_index(PhysiCell_constants::S_phase);
	int G2_index = flow_cytometry_separated_cycle_model.find_phase_index(PhysiCell_constants::G2_phase);
	int M_index = flow_cytometry_separated_cycle_model.find_phase_index(PhysiCell_constants::M_phase);

	// make this cell type randomly motile, less adhesive, greater survival, 
	// and less proliferative 
	
	epithelial = cell_defaults; 
	epithelial.type = 0; 
	epithelial.name = "epithelial cell"; 
	
	// make sure the new cell type has its own reference phenotype
	
	epithelial.parameters.pReference_live_phenotype = &( epithelial.phenotype ); 
	
	// enable random motility 
	epithelial.phenotype.motility.is_motile = true; 
	epithelial.phenotype.motility.persistence_time = parameters.doubles( "motile_cell_persistence_time" ); // 15.0; 
	epithelial.phenotype.motility.migration_speed = parameters.doubles( "motile_cell_migration_speed" ); // 0.25 micron/minute 
	epithelial.phenotype.motility.migration_bias = 0.0;// completely random 
	
	// Set cell-cell adhesion to 5% of other cells 
	epithelial.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "motile_cell_relative_adhesion" ); // 0.05; 
	
	// Set apoptosis to zero 
	epithelial.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "motile_cell_apoptosis_rate" ); // 0.0; 

	// Set proliferation to 10% of other cells. 
	// Alter the transition rate from G0G1 state to S state
	//motile_cell.phenotype.cycle.data.transition_rate(G0G1_index,S_index) *= 
	//	parameters.doubles( "motile_cell_relative_cycle_entry_rate" ); // 0.1; 
	
	//epithelial.phenotype.cycle.data.transition_rate(G0G1_index,S_index) *= 
	//	parameters.doubles(0.1); // 0.1;

	//epithelial.phenotype.cycle.data.transition_rate(M_index, G0G1_index) *= 
	//	parameters.doubles(0);

	///////////////////////////////////////////////////////////////////////////
	mesenchymal = cell_defaults; 
	mesenchymal.type = 1; 
	mesenchymal.name = "mesenchymal cell"; 
	
	// make sure the new cell type has its own reference phenotype
	
	mesenchymal.parameters.pReference_live_phenotype = &(mesenchymal.phenotype ); 
	
	// enable random motility 
	mesenchymal.phenotype.motility.is_motile = true; 
	mesenchymal.phenotype.motility.persistence_time = parameters.doubles( "motile_cell_persistence_time" ); // 15.0; 
	mesenchymal.phenotype.motility.migration_speed = parameters.doubles( "motile_cell_migration_speed" ); // 0.25 micron/minute 
	mesenchymal.phenotype.motility.migration_bias = 0.0;// completely random 
	
	// Set cell-cell adhesion to 5% of other cells 
	mesenchymal.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "motile_cell_relative_adhesion" ); // 0.05; 
	
	// Set apoptosis to zero 
	mesenchymal.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "motile_cell_apoptosis_rate" ); // 0.0; 

	// Set proliferation to 0% of other cells. 
	mesenchymal.phenotype.cycle.data.transition_rate(S_index, G2_index) *= 
		parameters.doubles(1);
	mesenchymal.phenotype.cycle.data.transition_rate(G2_index,M_index) *= 
		parameters.doubles(1);


	///////////////////////////////////////////////////////////////////////////
	stem_like = cell_defaults; 
	stem_like.type = 2; 
	stem_like.name = "stem_like cell"; 
	
	// make sure the new cell type has its own reference phenotype
	
	stem_like.parameters.pReference_live_phenotype = &( stem_like.phenotype ); 
	
	// enable random motility 
	stem_like.phenotype.motility.is_motile = true; 
	stem_like.phenotype.motility.persistence_time = parameters.doubles( "motile_cell_persistence_time" ); // 15.0; 
	stem_like.phenotype.motility.migration_speed = parameters.doubles( "motile_cell_migration_speed" ); // 0.25 micron/minute 
	stem_like.phenotype.motility.migration_bias = 0.0;// completely random 
	
	// Set cell-cell adhesion to 5% of other cells 
	stem_like.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "motile_cell_relative_adhesion" ); // 0.05; 
	
	// Set apoptosis to zero 
	stem_like.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "motile_cell_apoptosis_rate" ); // 0.0; 

	// Set proliferation to 0% of other cells. 
	stem_like.phenotype.cycle.data.transition_rate(S_index, G2_index) *= 
		parameters.doubles(1);
	stem_like.phenotype.cycle.data.transition_rate(G2_index,M_index) *= 
		parameters.doubles(1);


	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 

	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	/* now this is in XML 	
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	*/
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector(1 , 38.0); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 

	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 

	microenvironment.add_density( "glucose", "dimensionless" );
	microenvironment.diffusion_coefficients[1] = 1e3;
	microenvironment.decay_rates[1] = 0.05625;

	/* If we want to add drug treatment later
	microenvironment.add_density( "drug", "dimensionless" );
	microenvironment.diffusion_coefficients[2] = 1e3;
	microenvironment.decay_rates[2] = 0.15625;
	*/
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;

	for (int n = 0; n < 75; n++) {
		double x = rand() % 1800;
		double y = rand() % 1800;
		pC = create_cell(epithelial); 
		pC->assign_position( x-900, y-900, 0.0);

		x = rand() % 1800;
		y = rand() % 1800;
		pC = create_cell(mesenchymal);
		pC->assign_position( x-900, y-900, 0.0);

		x = rand() % 1800;
		y = rand() % 1800;
		pC = create_cell(stem_like);
		pC->assign_position( x-900, y-900, 0.0);
  	}
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell);
	
	//&& pCell->type == 0

	// Alive - Epithelial - Green
	if( pCell->phenotype.death.dead == false && pCell->type == 0)
	{
		 output[0] = "black"; 
		 output[2] = "rgb(0,180,0)"; 
	}

	// Alive - Mesenchymal - Blue
	if( pCell->phenotype.death.dead == false && pCell->type == 1)
	{
		 output[0] = "rgb(65,105,225)"; 
		 output[2] = "rgb(0,191,255)"; 
	}

	// Alive - Stem Like - Yellow
	if( pCell->phenotype.death.dead == false && pCell->type == 2)
	{
		 output[0] = "black"; 
		 output[2] = "rgb(180,180,0)"; 
	}


	// Apoptotic - Red
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )
	{
		output[0] = "rgb(125,0,0)";
		output[2] = "rgb(255,0,0)";
	}
	
	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "black";
		output[2] = "rgb(150,75,0)";
	}	

	return output; 
}
