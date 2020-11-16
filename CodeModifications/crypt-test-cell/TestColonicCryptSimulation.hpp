#ifndef TESTCOLONICCRYPTSIMULATION_HPP_
#define TESTCOLONICCRYPTSIMULATION_HPP_

// The first two header files are needed for parallelization.

#include <mpi.h>
#include <stdio.h>

#include <iostream>
#include <cxxtest/TestSuite.h>
#include <vector>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"

// The next header file defines a helper class for generating cells for crypt simulations. 

#include "CryptCellsGenerator.hpp"

/* The next header file defines a WntCellCycleModel, where the proliferative behaviour of a cell is
 * dependent on the concentration of Wnt at that point in space. Cells proliferate where there is a plentiful level of Wnt
 * and cease proliferation below a given threshold. */

#include "SimpleWntCellCycleModel.hpp"

/* The next header file defines a helper class for generating a suitable triangular mesh
 * for the crypt simulation, such that the cell corresponding to each node is initially
 * in mechanical equilibrium with its neighours and periodic boundary conditions are applied
 * at the left- and right-hand sides of the mesh (hence the "cylindrical"). */

#include "CylindricalHoneycombMeshGenerator.hpp"

/* The next header file defines a CellPopulation class that uses a triangular mesh, and allows
 * for the inclusion of 'ghost nodes'. These are nodes in the mesh that do not correspond
 * to cells, instead they help ensure that a sensible Delaunay triangulation is generated
 * at each timestep. This is because the triangulation algorithm requires a convex hull. */

#include "MeshBasedCellPopulationWithGhostNodes.hpp"

/* The next header file defines a force law, based on a linear spring, for describing
 * the mechanical interactions between neighbouring cells in the crypt. */

#include "GeneralisedLinearSpringForce.hpp"

/* The next header file defines the class that simulates the evolution of a CellPopulation,
 * specialized to deal with the cylindrical crypt geometry. */

#include "CryptSimulation2d.hpp"

/* The next header file defines a Wnt singleton class, which (if used) deals with the
 * imposed Wnt gradient in our crypt model. This affects cell proliferation in the case
 * where we construct each cell with a WntCellCycleModel. */

#include "WntConcentration.hpp"

/* The next header file defines a cell killer class, which implements sloughing of cells
 * into the lumen once they reach the top of the crypt. */

#include "SloughingCellKiller.hpp"

// The next header file defines another cell killer class, which implements the homeostatic apoptosis of cells.

#include "HomeostaticCellKiller.hpp"

// The next header file defines another cell killer class, which simulates crypt fission.

#include "CryptFissionCellKiller.hpp"

// These headers are used for defining and recording mutations.

#include "WildTypeCellMutationState.hpp"
#include "MMRTwoHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "BetaCateninTwoHitCellMutationState.hpp"
#include "MMR_ApcOneHitCellMutationState.hpp"
#include "MMR_ApcTwoHitCellMutationState.hpp"
#include "MMR_BCOneHitCellMutationState.hpp"
#include "MMR_BCTwoHitCellMutationState.hpp"
#include "MMR_BCOneHit_ApcOneHitCellMutationState.hpp"
#include "BCOneHit_ApcOneHitCellMutationState.hpp"
#include "ApcLOHCellMutationState.hpp"
#include "BCOneHit_ApcLOHCellMutationState.hpp"
#include "MMR_ApcLOHCellMutationState.hpp"
#include "MMR_BCOneHit_ApcLOHCellMutationState.hpp"
#include "BCLOHApcOneHitCellMutationState.hpp"
#include "BCLOHMMRTwoHitCellMutationState.hpp"
#include "BCLOHMMR_ApcOneHitCellMutationState.hpp"
#include "BCLOHApcLOHCellMutationState.hpp"
#include "BCLOHMMR_ApcLOHCellMutationState.hpp"
#include "BCLOHCellMutationState.hpp"

// The next header allows us to track how many cells of each mutation state there are.

#include "CellMutationStatesCountWriter.hpp"

// The final header tracks the volume of all cells, which is needed for apoptosis.

#include "VolumeTrackingModifier.hpp"


// Next, we define the test class, which inherits from AbstractCellBasedTestSuite.

class TestColonicCryptSimulation : public AbstractCellBasedTestSuite
{
public:
  	

    void TestColonicCrypt()
    {
		/* We create boost shared pointers to any mutations we wish to use.
	     * We need to do this using the CellPropertyRegistry, otherwise
	     * the numbers of each type of mutation aren't correctly tracked.
		 * Here, we include all (relevant) combinations of MMR, CTNNB1/beta-catenin (BC), and APC
		 * (point mutation and LOH) mutations, and a wild type (WT) mutation state.
	     * Each mutation has a different effect on the cell cycle model, see SimpleWntCellCycleModel.cpp
		 * for details. */

		boost::shared_ptr<AbstractCellProperty> p_state_0(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
	    boost::shared_ptr<AbstractCellProperty> p_state_1(CellPropertyRegistry::Instance()->Get<MMRTwoHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_2(CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_3(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_4(CellPropertyRegistry::Instance()->Get<BetaCateninOneHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_5(CellPropertyRegistry::Instance()->Get<MMRApcOneHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_6(CellPropertyRegistry::Instance()->Get<MMRApcTwoHitCellMutationState>());		
		boost::shared_ptr<AbstractCellProperty> p_state_7(CellPropertyRegistry::Instance()->Get<MMRBCOneHitCellMutationState>());	
		boost::shared_ptr<AbstractCellProperty> p_state_8(CellPropertyRegistry::Instance()->Get<BCOneHitApcOneHitCellMutationState>());	
		boost::shared_ptr<AbstractCellProperty> p_state_9(CellPropertyRegistry::Instance()->Get<MMRBCOneHitApcOneHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_10(CellPropertyRegistry::Instance()->Get<BetaCateninTwoHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_11(CellPropertyRegistry::Instance()->Get<MMRBCTwoHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_12(CellPropertyRegistry::Instance()->Get<ApcLOHCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_13(CellPropertyRegistry::Instance()->Get<MMRApcLOHCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_14(CellPropertyRegistry::Instance()->Get<BCOneHitApcLOHCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_15(CellPropertyRegistry::Instance()->Get<MMRBCOneHitApcLOHCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_16(CellPropertyRegistry::Instance()->Get<BCLOHMMRTwoHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_17(CellPropertyRegistry::Instance()->Get<BCLOHApcOneHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_18(CellPropertyRegistry::Instance()->Get<BCLOHMMRApcOneHitCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_19(CellPropertyRegistry::Instance()->Get<BCLOHApcLOHCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_20(CellPropertyRegistry::Instance()->Get<BCLOHMMRApcLOHCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_state_21(CellPropertyRegistry::Instance()->Get<BCLOHCellMutationState>());
		
		/* We need to distinguish between MLH1 and MSH2 germline mutation carriers, since MLH1 LOH of the WT allele can induce a CTNNB1 double hit,
		 * or lead to the loss of a CTNNB1 mutation.
		 * To simulate different patients, we need to alter the boolean variable here and in SimpleWntCellCycleModel.cpp. */

		bool MLH1 = true; // set germline mutation: true = MLH1, false = MSH2.

		/* We can sequentially simulate crypts by using a for-loop.
		 * We define the number of initially simulated crypts, given by the last integer of the for-loop.
		 * We also initialize the total number of crypts, later calculated by adding all crypts from the simulation.
		 * At the end, we want to know their difference, i.e. the number of total crypt fissions.
		 * We further define the number of MMR-deficient crypts, which interests us regarding the data.
		 * We can trace back the distribution to the foci, whose number we also calculate, via the output. */

		int initial_crypt_no = 1;
		int total_crypt_no = 0;
		int MMR_crypt_no = 0;
		int MMR_DCF_no = 0;
		
		// Initialize the MPI environment, to simulate multiple crypts simultaneously. We collectively call the crypts simultated on each processor "clusters".

		MPI_Init(NULL, NULL);
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		
		
		for (int j=1; j<2; j++)
		{

		    /* We set the default height of the crypt. As well as passing this variable into the SloughingCellKiller,
		     * we will pass it to the WntConcentration object (see below).
			 * Further we set the baseline probability of crypt fission withing an hour, which is passed into the CryptFissionCellKiller,
			 * where it is increased for MMR-deficient crypts. */
		    
		    double crypt_height = 69.0;
			double mCryptFissionRate = 0.0000012627;  // per 1 hour
			
			/* We initialize the crypt counter, which is incremented whenever crypt fission occurs. The crypt fission rate is multiplied accordingly.
			 * This variable is also passed into the CryptFissionCellKiller. */
			 
			int crypt_counter = 1;

		    // We first create a cylindrical mesh, and get the cell location indices.
		   
		    CylindricalHoneycombMeshGenerator generator(20, 80, 2); // (circumference, height, ghost node layers)
		    Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

		    std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		    // We create the cells. Here we use a SimpleWntCellCycleModel.

		    std::vector<CellPtr> cells;
		    CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
		    cells_generator.Generate(cells, p_mesh, location_indices, true);

			// We create the cell population.

		    MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices, 14.0);

		    // In order to visualize mutant cells and to count how many cells there are of each type we need to use the following command.

		    cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
		

		    /* When using a SimpleWntCellCycleModel, we need a way of telling each cell what the Wnt concentration
		     * is at its location. To do this, we set up a WntConcentration object. Like SimulationTime,
		     * WntConcentration is a singleton class, so when instantiated it is accessible from anywhere in
		     * the code (and in particular, all cells and cell-cycle models can access it). We need to say what
		     * the profile of the Wnt concentation should be up the crypt: here, we say it is LINEAR (linear
		     * decreasing from 1 to 0 from the bottom of the crypt to the top). We also need to inform the
		     * WntConcentration of the cell population and the height of the crypt. */

		    WntConcentration<2>::Instance()->SetType(LINEAR);
		    WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
		    WntConcentration<2>::Instance()->SetCryptLength(crypt_height);

			/* We want to obtain different outcomes in each cluster. Therefore we reseed the RandomNumberGenerator to a number 
			 * which is very likely to be new. Here we use getpid(), which returns the system's process ID for the current program. */

			RandomNumberGenerator::Instance()->Reseed(getpid()); 
			 
		    // Create a simulator.

		    CryptSimulation2d simulator(cell_population);

			// Create a specific output directory for every crypt, to be able to view specific simulations and results.
	
			std::string outputdirectory;
			outputdirectory = "Crypt" + std::to_string(rank) + "." + std::to_string(j);
		    simulator.SetOutputDirectory(outputdirectory);
		    simulator.SetSamplingTimestepMultiple(450); // only generate output every x*(1/45) hours.

		    // We create a force law, the three cell killers and the volume tracker, and pass these objects to the simulator.

		    MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
		    simulator.AddForce(p_linear_force);
		    MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
		    simulator.AddCellKiller(p_killer);
			MAKE_PTR_ARGS(HomeostaticCellKiller<2>, p_killer1, (&cell_population));
		    simulator.AddCellKiller(p_killer1);
			MAKE_PTR_ARGS(CryptFissionCellKiller<2>, p_killer2, (&cell_population, mCryptFissionRate, crypt_counter));
		    simulator.AddCellKiller(p_killer2);
			MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        	simulator.AddSimulationModifier(p_modifier);
			
			
			// Set up boolean variables to track mutation states of the stem cell.

			bool Is_SC_MMR_mut = false;
			bool Is_SC_MMR_LOH_mut = false; // one LOH event for MLH1, two LOH events for MSH2
			bool Is_SC_MLH1_LOH_mut = false; // two MLH1 LOH events
			bool Is_SC_APC_One_Hit_LOH_mut = false;
			bool Is_SC_APC_One_Hit_pt_mut = false;
			bool Is_SC_APC_Two_Hit_mut = false;
			bool Is_SC_BC_One_Hit_mut = false;
			bool Is_SC_BC_Two_Hit_mut = false;
			bool Is_BC_affected = false; // LOH of CTNNB1. Does not need a proper mutation state, because no function is changed.

			int SC_number = 5; // number of stem cells in the crypt -1 (e.g. 5 for 6 stem cells)
			int SC_index = 0; // number of stem cells that have once populated the crypt, ranging from 0 to SC_number
			int current_SC = 0; // number of the stem cell currently populating the crypt, also ranging from 0 to SC_number
			std::vector<std::vector<bool> > all_mut_states(SC_number+1); // 2D-vector that stores the mutation combination of each stem cell
			std::vector<bool> current_mut_state(9,false); // vector that stores the mutation combination of the current stem cell

			// APC-mutated stem cell as initial condition.

			Is_SC_APC_One_Hit_LOH_mut = true;
			for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
		         	cell_iter != cell_population.End();
		         	++cell_iter)
		    {
		        double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

		       	if (cell_height <= 0.5) // all cells of the lowest TA cell row are mutated. The new mutation status depends on previous mutations.
		       	{
						cell_iter->SetMutationState(p_state_12); // APC LOH mutation
				}
			}

			// Output that we start the simulation.

			cout << "Starting simulation of crypt " << j << " in cluster " << rank << endl;

			// We simulate in multiples of the stem cell cycle time. The integer i refers to the number of cycles.

			for (int i=1; i<7; i++)
			{

				simulator.SetEndTime(i*1680); // The number which is multiplied by i is the SC cycle time (in hours).
				simulator.Solve();
		
				// End of the stem cell cycle
				
				cout << "End of SC cycle " << i << " in crypt " << j << " in cluster " << rank << "." << endl;
				
				// We calculate the percentage of MMR-deficient cells within the crypt.

				double mmr_percentage = (double)cell_population.GetNumMMRCells()/(double)cell_population.GetNumRealCells();	

				if (mmr_percentage > 0.8) // declare crypt MMR-deficient if over 80% of its cells are
				{
					cout << "Crypt " << j << " in cluster " << rank << " is MMR-deficient after " << i << " SC cycles." << endl;  // output if crypt is MMR-deficient, just as information
				}

				current_mut_state[0]=Is_SC_MMR_mut;
				current_mut_state[1]=Is_SC_MMR_LOH_mut;
				current_mut_state[2]=Is_SC_MLH1_LOH_mut;
				current_mut_state[3]=Is_SC_APC_One_Hit_LOH_mut;
				current_mut_state[4]=Is_SC_APC_One_Hit_pt_mut;
				current_mut_state[5]=Is_SC_APC_Two_Hit_mut;
				current_mut_state[6]=Is_SC_BC_One_Hit_mut;
				current_mut_state[7]=Is_SC_BC_Two_Hit_mut;
				current_mut_state[8]=Is_BC_affected;
				all_mut_states[current_SC] = current_mut_state; // store the mutation combination of each stem cell

				// We first consider the loss of stem cells due to lethal mutations occurring upon cell division.

				double SC_loss_prob = 0.0001; // Probability that a stem cell dies (corresponds to about 1/30 % deleterious mutations)
				
				if (Is_SC_MMR_mut)
				{									
					SC_loss_prob *= 100.0; // Increase SC death probability 100-fold
				}

				//SC_loss_prob *= 100.0; // mice simulation


				if (RandomNumberGenerator::Instance()->ranf() < SC_loss_prob) // stem cell is lost
				{	

					cout << "SC " << current_SC+1 << " of crypt " << j << " in cluster " << rank  <<  " is lost at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of SC loss
					
					// The stem cell loss is compensated for by symmetric division of an adjacent stem cell.				
					
					if (current_SC < SC_number)
					{
						all_mut_states[current_SC] = all_mut_states[current_SC+1];
						cout << "The loss is replenished by SC " << current_SC+2 << "." << endl;
					}
					else
					{
						all_mut_states[current_SC] = all_mut_states[0];
						cout << "The loss is replenished by SC " << 1 << "." << endl;
					}

					// A new stem cell is chosen for populating the crypt in the meantime.

					// Change cell numbering to make things easier, lost cell is now the "last" (i.e. the cell index) cell (in case it was not already).
					std::vector<bool> temp = all_mut_states[SC_index]; // store status of the current "last" cell in the ordering temporarily
					all_mut_states[SC_index] = all_mut_states[current_SC]; // first swap
					all_mut_states[current_SC]= temp; // second swap
				
					if (SC_index != current_SC)
					{
						cout << "We changed the cell numbering. SC " << current_SC+1 << " is now SC " << SC_index+1 << " and vice versa." << endl;
					}

					double p_old_SC = (double)SC_index/(double)SC_number; // the probability that the next stem cell is an "old" one

					if (RandomNumberGenerator::Instance()->ranf() < p_old_SC) // old stem cell
					{
						for (int i=1; i < SC_index; i++)
						{
							// replace mutation status by the ones of the old stem cell						
	
							if (RandomNumberGenerator::Instance()->ranf() < (double)(SC_index - i)/(double)(SC_index +1 -i))
							{
								Is_SC_MMR_mut = all_mut_states[i][0];
								Is_SC_MMR_LOH_mut = all_mut_states[i][1];
								Is_SC_MLH1_LOH_mut = all_mut_states[i][2];
								Is_SC_APC_One_Hit_LOH_mut = all_mut_states[i][3];
								Is_SC_APC_One_Hit_pt_mut = all_mut_states[i][4];
								Is_SC_APC_Two_Hit_mut = all_mut_states[i][5];
								Is_SC_BC_One_Hit_mut = all_mut_states[i][6];
								Is_SC_BC_Two_Hit_mut = all_mut_states[i][7];
								Is_BC_affected = all_mut_states[i][8];
								if (i == (SC_index-1)) // only last one is relevant, i.e. the "previous" cell
								{
									cout << "It is now populated by SC " << i+1 << " again." << endl;
									current_SC = i;
								}
							}
							else
							{
								Is_SC_MMR_mut = all_mut_states[i-1][0];
								Is_SC_MMR_LOH_mut = all_mut_states[i-1][1];
								Is_SC_MLH1_LOH_mut = all_mut_states[i-1][2];
								Is_SC_APC_One_Hit_LOH_mut = all_mut_states[i-1][3];
								Is_SC_APC_One_Hit_pt_mut = all_mut_states[i-1][4];
								Is_SC_APC_Two_Hit_mut = all_mut_states[i-1][5];
								Is_SC_BC_One_Hit_mut = all_mut_states[i-1][6];
								Is_SC_BC_Two_Hit_mut = all_mut_states[i-1][7];
								Is_BC_affected = all_mut_states[i-1][8];
								cout << "It is now populated by SC " << i << " again." << endl;
								current_SC = i-1;
								break;
							}
						}
					}
					else // new stem cell
					{
						SC_index++;
						current_SC = SC_index;
						cout << "It is now populated by a new SC, namely SC " << current_SC+1 << ". Mutations are lost." << endl; 
					
						//Reset boolean variables
						Is_SC_MMR_mut = false;
						Is_SC_MMR_LOH_mut = false;
						Is_SC_MLH1_LOH_mut = false;
						Is_SC_APC_One_Hit_LOH_mut = false;
						Is_SC_APC_One_Hit_pt_mut = false;
						Is_SC_APC_Two_Hit_mut = false;
						Is_SC_BC_One_Hit_mut = false;
						Is_SC_BC_Two_Hit_mut = false;
						Is_BC_affected = false;

						for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				     	cell_iter != cell_population.End();
				     	++cell_iter)
						{
				    		double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

				    		if (cell_height <= 0.5) 
				    		{
								cell_iter->SetMutationState(p_state_0); // all cells of the lowest TA cell row lose their mutations.
				    		}
						}
					}
				}
				
				// Now, consider a change of the stem cell which populates the crypt.

				double SC_change_prob = 0.5; // Probability that a different SC populates the crypt (should occur every 5 months on average)

				if (RandomNumberGenerator::Instance()->ranf() < SC_change_prob) 
				{	
					// Change cell numbering to make things easier, changed cell is now the "last" (i.e. the cell index) cell (in case it was not already).
					std::vector<bool> temp = all_mut_states[SC_index]; // store status of the "last" cell in the ordering temporarily
					all_mut_states[SC_index] = all_mut_states[current_SC]; // first swap
					all_mut_states[current_SC]= temp; // second swap

					if (SC_index != current_SC)
					{
						cout << "We changed the cell numbering. SC " << current_SC+1 << " is now SC " << SC_index+1 << " and vice versa." << endl;
					}
					
					double p_old_SC = (double)SC_index/(double)SC_number; // the probability that the next stem cell is an "old" one

					if (RandomNumberGenerator::Instance()->ranf() < p_old_SC) // old stem cell
					{
						for (int i=1; i < SC_index; i++)
						{
							// replace mutation status by the ones of the old stem cell

							if (RandomNumberGenerator::Instance()->ranf() < (double)(SC_index - i)/(double)(SC_index +1 -i))
							{
								Is_SC_MMR_mut = all_mut_states[i][0];
								Is_SC_MMR_LOH_mut = all_mut_states[i][1];
								Is_SC_MLH1_LOH_mut = all_mut_states[i][2];
								Is_SC_APC_One_Hit_LOH_mut = all_mut_states[i][3];
								Is_SC_APC_One_Hit_pt_mut = all_mut_states[i][4];
								Is_SC_APC_Two_Hit_mut = all_mut_states[i][5];
								Is_SC_BC_One_Hit_mut = all_mut_states[i][6];
								Is_SC_BC_Two_Hit_mut = all_mut_states[i][7];
								Is_BC_affected = all_mut_states[i][8];
								if (i == (SC_index-1)) // only last one is relevant, i.e. the "previous" cell
								{
									cout << "Stem Cell Change: Crypt " << j << " in cluster " << rank << " is now populated by SC " << i+1 << " again." << endl;
									current_SC = i;
								}
							}
							else
							{
								Is_SC_MMR_mut = all_mut_states[i-1][0];
								Is_SC_MMR_LOH_mut = all_mut_states[i-1][1];
								Is_SC_MLH1_LOH_mut = all_mut_states[i-1][2];
								Is_SC_APC_One_Hit_LOH_mut = all_mut_states[i-1][3];
								Is_SC_APC_One_Hit_pt_mut = all_mut_states[i-1][4];
								Is_SC_APC_Two_Hit_mut = all_mut_states[i-1][5];
								Is_SC_BC_One_Hit_mut = all_mut_states[i-1][6];
								Is_SC_BC_Two_Hit_mut = all_mut_states[i-1][7];
								Is_BC_affected = all_mut_states[i-1][8];
								cout << "Stem Cell Change: Crypt " << j << " in cluster " << rank << " is now populated by SC " << i << " again." << endl;
								current_SC = i-1;
								break;
							}
						}
					}
					else // new stem cell
					{
						SC_index++;
						current_SC = SC_index;
						cout << "Crypt " << j << " in cluster " << rank  << " is populated by a new SC, namely SC " << current_SC+1 << ". Mutations are lost." << endl;
					
						// Reset boolean variables
						Is_SC_MMR_mut = false;
						Is_SC_MMR_LOH_mut = false;
						Is_SC_MLH1_LOH_mut = false;
						Is_SC_APC_One_Hit_LOH_mut = false;
						Is_SC_APC_One_Hit_pt_mut = false;
						Is_SC_APC_Two_Hit_mut = false;
						Is_SC_BC_One_Hit_mut = false;
						Is_SC_BC_Two_Hit_mut = false;
						Is_BC_affected = false;

						for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				     	cell_iter != cell_population.End();
				     	++cell_iter)
						{
				    		double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

				    		if (cell_height <= 0.5) 
				    		{
								cell_iter->SetMutationState(p_state_0); // all cells of the lowest TA cell row lose their mutations.
				    		}
						}
					}
				}
				
				/* After each stem cell cycle, the SC can mutate according to certain probabilities.
				 * We simulate a SC mutation by mutating all cells from the lowest row of the crypt. */

				// We first define the mutation probabilities. 

				double MMRMutProb; // Probability of a sporadic MMR mutation. Depends on the germline mutation of the patient.

				if (MLH1)
				{
					MMRMutProb = 0.000003546875;
				}
				else
				{
					MMRMutProb = 0.000004375; // MSH2
				}
			
				double MMRLOHProb;

				if (MLH1)
				{
					MMRLOHProb = 4.0*MMRMutProb; // consider both alleles for MLH1, due to importance for CTNNB1 LOH
				}
				else
				{
					MMRLOHProb = 0.00000987; // MSH2, only consider second allele
				}
		
				double ApcHitProb = 0.0000075; // Probability of a sporadic Apc mutation.
				double ApcLOHProb = 0.0000343;
				double BetaCateninProb = 0.000000003125; // Probability of a sporadic Beta-Catenin mutation.
				double BetaCateninLOHProb = 0.000010116; // Probability of Beta-Catenin LOH (important if second allele becomes mutated)

				if (Is_SC_MMR_mut)
				{
					// Increase the probability of further point mutations 100-fold.
					ApcHitProb *= 100.0;
					BetaCateninProb *= 100.0;
					MMRMutProb = 0.0;
				}

				/*// mice simulation: all cells are MMR-deficient, which here is the wild-type
				ApcHitProb *= 100.0;
				BetaCateninProb *= 100.0;
				MMRMutProb = 0.0;
				MMRLOHProb = 0.0;*/

				if (Is_SC_MMR_LOH_mut)
				{
					if (MLH1)
					{
						MMRLOHProb *= 0.5;
					}
					else
					{
						MMRLOHProb = 0.0; // do not care about further LOH events in the case of MSH2
					}
				}

				if (Is_SC_MLH1_LOH_mut)
				{
					MMRLOHProb = 0.0;
				}

				if (Is_SC_APC_One_Hit_LOH_mut)
				{
					ApcLOHProb *= 0.5; 
					ApcHitProb *= 0.5;
				}

				if (Is_SC_APC_One_Hit_pt_mut)
				{
					ApcHitProb *= 0.5;
				}

				if (Is_SC_APC_Two_Hit_mut)
				{
					ApcHitProb = 0.0;
				}

				if (Is_SC_BC_One_Hit_mut)
				{
					BetaCateninProb *= 0.5;
				}

				if (Is_SC_BC_Two_Hit_mut)
				{
					BetaCateninProb = 0.0;
				}

				if (Is_BC_affected)
				{
					BetaCateninProb *= 0.5; // There now is one less WT allele to be mutated.
					BetaCateninLOHProb *= 0.5;
				}


				// We start with MMR mutations. LOH is first, which is the more complicated case.

				if (RandomNumberGenerator::Instance()->ranf() < MMRLOHProb)
				{	
					cout << "SC " << current_SC+1 << " undergoes MMR LOH in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output time of event

					bool NewStemCell = false; // track whether the stem cell dies. If so, then the MLH1 mutation does not matter anymore.

					// First, take care of concurrent CTNNB1 loss for MLH1

					if (MLH1 && RandomNumberGenerator::Instance()->ranf() < 0.8) // MLH1 LOH affects CTNNB1 (in 80% of cases)
					{
									
						if (Is_BC_affected && RandomNumberGenerator::Instance()->ranf() < 0.5) // one allele was already lost and second one now too
						{	
							NewStemCell = true;
							// stem cell dies, copy from above:

							cout << "SC " << current_SC+1 << " of crypt " << j << " in cluster " << rank  <<  " is lost due to CTNNB1 LOH at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of SC loss
					
							// The stem cell loss is compensated for by symmetric division of an adjacent stem cell.				
							
							if (current_SC < SC_number)
							{
								all_mut_states[current_SC] = all_mut_states[current_SC+1];
								cout << "The loss is replenished by SC " << current_SC+2 << "." << endl;
							}
							else
							{
								all_mut_states[current_SC] = all_mut_states[0];
								cout << "The loss is replenished by SC " << 1 << "." << endl;
							}

							// A new stem cell is chosen for populating the crypt in the meantime.

							// Change cell numbering to make things easier, lost cell is now the "last" (i.e. the cell index) cell (in case it was not already).
							std::vector<bool> temp = all_mut_states[SC_index]; // store status of the current "last" cell in the ordering temporarily
							all_mut_states[SC_index] = all_mut_states[current_SC]; // first swap
							all_mut_states[current_SC]= temp; // second swap
						
							if (SC_index != current_SC)
							{
								cout << "We changed the cell numbering. SC " << current_SC+1 << " is now SC " << SC_index+1 << " and vice versa." << endl;
							}

							double p_old_SC = (double)SC_index/(double)SC_number; // the probability that the next stem cell is an "old" one

							if (RandomNumberGenerator::Instance()->ranf() < p_old_SC) // old stem cell
							{
								for (int i=1; i < SC_index; i++)
								{
									// replace mutation status by the ones of the old stem cell						
			
									if (RandomNumberGenerator::Instance()->ranf() < (double)(SC_index - i)/(double)(SC_index +1 -i))
									{
										Is_SC_MMR_mut = all_mut_states[i][0];
										Is_SC_MMR_LOH_mut = all_mut_states[i][1];
										Is_SC_MLH1_LOH_mut = all_mut_states[i][2];
										Is_SC_APC_One_Hit_LOH_mut = all_mut_states[i][3];
										Is_SC_APC_One_Hit_pt_mut = all_mut_states[i][4];
										Is_SC_APC_Two_Hit_mut = all_mut_states[i][5];
										Is_SC_BC_One_Hit_mut = all_mut_states[i][6];
										Is_SC_BC_Two_Hit_mut = all_mut_states[i][7];
										Is_BC_affected = all_mut_states[i][8];
										if (i == (SC_index-1)) // only last one is relevant, i.e. the "previous" cell
										{
											cout << "It is now populated by SC " << i+1 << " again." << endl;
											current_SC = i;
										}
									}
									else
									{
										Is_SC_MMR_mut = all_mut_states[i-1][0];
										Is_SC_MMR_LOH_mut = all_mut_states[i-1][1];
										Is_SC_MLH1_LOH_mut = all_mut_states[i-1][2];
										Is_SC_APC_One_Hit_LOH_mut = all_mut_states[i-1][3];
										Is_SC_APC_One_Hit_pt_mut = all_mut_states[i-1][4];
										Is_SC_APC_Two_Hit_mut = all_mut_states[i-1][5];
										Is_SC_BC_One_Hit_mut = all_mut_states[i-1][6];
										Is_SC_BC_Two_Hit_mut = all_mut_states[i-1][7];
										Is_BC_affected = all_mut_states[i-1][8];
										cout << "It is now populated by SC " << i << " again." << endl;
										current_SC = i-1;
										break;
									}
								}
							}
							else // new stem cell
							{
								SC_index++;
								current_SC = SC_index;
								cout << "It is now populated by a new SC, namely SC " << current_SC+1 << ". Mutations are lost." << endl; 
							
								//Reset boolean variables
								Is_SC_MMR_mut = false;
								Is_SC_MMR_LOH_mut = false;
								Is_SC_MLH1_LOH_mut = false;
								Is_SC_APC_One_Hit_LOH_mut = false;
								Is_SC_APC_One_Hit_pt_mut = false;
								Is_SC_APC_Two_Hit_mut = false;
								Is_SC_BC_One_Hit_mut = false;
								Is_SC_BC_Two_Hit_mut = false;
								Is_BC_affected = false;

								for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
							 	cell_iter != cell_population.End();
							 	++cell_iter)
								{
									double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

									if (cell_height <= 0.5) 
									{
										cell_iter->SetMutationState(p_state_0); // all cells of the lowest TA cell row lose their mutations.
									}
								}
							}
							// end of copy
						}
						else if (!Is_BC_affected) // no CTNNB1 LOH so far
						{
							Is_BC_affected = true; // BC is lost
							BetaCateninProb *= 0.5; // There now is one less WT allele to be mutated.
							BetaCateninLOHProb *= 0.5;

							if (Is_SC_BC_One_Hit_mut) // one point mut had already occurred
							{
								if (RandomNumberGenerator::Instance()->ranf() < 0.5) // BC mutation is on other allele, now we have a biallelic mutation
								{
									cout << "SC " << current_SC+1 << " undergoes beta-catenin double hit due to MLH1 LOH in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
									Is_SC_BC_One_Hit_mut = false;
									Is_SC_BC_Two_Hit_mut = true;
									BetaCateninProb = 0.0;
								}
								else // BC mutation is on lost allele, now we lose the point mutation
								{	
									cout << "Beta-catenin mutation is lost in SC " << current_SC+1 << " due to MLH1 LOH in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
									BetaCateninProb *= 2.0; // correct factor from above, since probability was already halved
									Is_SC_BC_One_Hit_mut = false; // mutation is lost
								}
							}
						}
					}
	

					// Now take care of the MMR LOH itself, if it is still the same stem cell.

					if (!NewStemCell)
					{
						if (Is_SC_MMR_LOH_mut) // second LOH event, must be MLH1, otherwise LOH probability would have been 0
						{
							Is_SC_MLH1_LOH_mut = true; // no further LOH events can occur
						}
						else // first LOH event
						{
							Is_SC_MMR_LOH_mut = true;
						}

						if ((MLH1 && Is_SC_MMR_LOH_mut && !Is_SC_MLH1_LOH_mut && RandomNumberGenerator::Instance()->ranf() < 0.5) || (MLH1 && Is_SC_MLH1_LOH_mut) || !MLH1) // first LOH and WT allele, second LOH, MSH2 LOH
						{	
							Is_SC_MMR_mut = true; // SC is now MMR-deficient
			
							// Now all cells of the lowest TA cell row are assigned an MMR mutation. The new mutation state depends on previous mutations.

							for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
						 	cell_iter != cell_population.End();
						 	++cell_iter)
							{
								double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

								if (cell_height <= 0.5) // lowest row
								{
									if (cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_1); // MMR mutation
									}
									else if (cell_iter->GetMutationState()->IsType<BCLOHCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_16); // CTNNB1 LOH and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<ApcOneHitCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_5); // APC+- and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<ApcLOHCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_13); 
									}
									else if (cell_iter->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_18); // APC+-, CTNNB1 LOH and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_20);
									}
									else if (cell_iter->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_6); // APC-- and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
									{	
											cell_iter->SetMutationState(p_state_7); // beta-catenin + and MMR	
									}
									else if (cell_iter->GetMutationState()->IsType<BetaCateninTwoHitCellMutationState>())
									{					
										cell_iter->SetMutationState(p_state_11); // beta-catenin ++ and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>())
									{
											cell_iter->SetMutationState(p_state_9); // beta-catenin +, APC +- and MMR	
									}
									else if (cell_iter->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
									{
											cell_iter->SetMutationState(p_state_15); // beta-catenin +, APC +- and MMR	
									}
								}
							}
						}
					}
				}
									
				// Now do point mutations. If the stem cell already is MMR-deficient, we do not give any output and only mutate the cells again.

				if (RandomNumberGenerator::Instance()->ranf() < MMRMutProb || Is_SC_MMR_mut ) 
				{	
					if (!Is_SC_MMR_mut)	// relevant MMR point mutation occurs
					{
						cout << "SC " << current_SC+1 << " undergoes MMR point mutation in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
						Is_SC_MMR_mut = true; // cell is now MMR-deficient
					}
			
					for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
		         	cell_iter != cell_population.End();
		         	++cell_iter)
		    		{
		        		double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

		        		if (cell_height <= 0.5) // all cells of the lowest TA cell row are mutated. The new mutation status depends on previous mutations.
		        		{
							if (cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_1); // MMR mutation
							}
							else if (cell_iter->GetMutationState()->IsType<BCLOHCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_16); // CTNNB1 LOH and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<ApcOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_5); // APC+- and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<ApcLOHCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_13); 
							}
							else if (cell_iter->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_18); // APC+-, CTNNB1 LOH and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_20);
							}
							else if (cell_iter->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_6); // APC-- and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_7); // beta-catenin + and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<BetaCateninTwoHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_11); // beta-catenin ++ and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_9); // beta-catenin +, APC+- and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_15); 
							}
		        		}
					}
				}

				// Now do APC mutations. LOH again is first, only now we do not always want to know, since no other gene is involved.
				// We have to do it slightly differently than as done in the MMR case. First we compute whether a mutation occurs.
		
				if (RandomNumberGenerator::Instance()->ranf() < ApcLOHProb)
				{	

					bool NewStemCell = false; // again track whether stem cell dies			

					
					if (Is_SC_APC_One_Hit_LOH_mut) // second LOH event on the other allele (probability was halved)
					{
						NewStemCell = true; // cell dies

						// copy from above again:
						cout << "SC " << current_SC+1 << " of crypt " << j << " in cluster " << rank  <<  " is lost due to APC LOH at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of SC loss
					
						// The stem cell loss is compensated for by symmetric division of an adjacent stem cell.				
							
						if (current_SC < SC_number)
						{
							all_mut_states[current_SC] = all_mut_states[current_SC+1];
							cout << "The loss is replenished by SC " << current_SC+2 << "." << endl;
						}
						else
						{
							all_mut_states[current_SC] = all_mut_states[0];
							cout << "The loss is replenished by SC " << 1 << "." << endl;
						}

						// A new stem cell is chosen for populating the crypt in the meantime.

						// Change cell numbering to make things easier, lost cell is now the "last" (i.e. the cell index) cell (in case it was not already).
						std::vector<bool> temp = all_mut_states[SC_index]; // store status of the current "last" cell in the ordering temporarily
						all_mut_states[SC_index] = all_mut_states[current_SC]; // first swap
						all_mut_states[current_SC]= temp; // second swap
						
						if (SC_index != current_SC)
						{
							cout << "We changed the cell numbering. SC " << current_SC+1 << " is now SC " << SC_index+1 << " and vice versa." << endl;
						}

						double p_old_SC = (double)SC_index/(double)SC_number; // the probability that the next stem cell is an "old" one

						if (RandomNumberGenerator::Instance()->ranf() < p_old_SC) // old stem cell
						{
							for (int i=1; i < SC_index; i++)
							{
								// replace mutation status by the ones of the old stem cell						
			
								if (RandomNumberGenerator::Instance()->ranf() < (double)(SC_index - i)/(double)(SC_index +1 -i))
								{
									Is_SC_MMR_mut = all_mut_states[i][0];
									Is_SC_MMR_LOH_mut = all_mut_states[i][1];
									Is_SC_MLH1_LOH_mut = all_mut_states[i][2];
									Is_SC_APC_One_Hit_LOH_mut = all_mut_states[i][3];
									Is_SC_APC_One_Hit_pt_mut = all_mut_states[i][4];
									Is_SC_APC_Two_Hit_mut = all_mut_states[i][5];
									Is_SC_BC_One_Hit_mut = all_mut_states[i][6];
									Is_SC_BC_Two_Hit_mut = all_mut_states[i][7];
									Is_BC_affected = all_mut_states[i][8];
									if (i == (SC_index-1)) // only last one is relevant, i.e. the "previous" cell
									{
										cout << "It is now populated by SC " << i+1 << " again." << endl;
										current_SC = i;										
									}
								}
								else
								{
									Is_SC_MMR_mut = all_mut_states[i-1][0];
									Is_SC_MMR_LOH_mut = all_mut_states[i-1][1];
									Is_SC_MLH1_LOH_mut = all_mut_states[i-1][2];
									Is_SC_APC_One_Hit_LOH_mut = all_mut_states[i-1][3];
									Is_SC_APC_One_Hit_pt_mut = all_mut_states[i-1][4];
									Is_SC_APC_Two_Hit_mut = all_mut_states[i-1][5];
									Is_SC_BC_One_Hit_mut = all_mut_states[i-1][6];
									Is_SC_BC_Two_Hit_mut = all_mut_states[i-1][7];
									Is_BC_affected = all_mut_states[i-1][8];
									cout << "It is now populated by SC " << i << " again." << endl;
									current_SC = i-1;
									break;
								}
							}
						}
						else // new stem cell
						{
							SC_index++;
							current_SC = SC_index;
							cout << "It is now populated by a new SC, namely SC " << current_SC+1 << ". Mutations are lost." << endl; 
						
							//Reset boolean variables
							Is_SC_MMR_mut = false;
							Is_SC_MMR_LOH_mut = false;
							Is_SC_MLH1_LOH_mut = false;
							Is_SC_APC_One_Hit_LOH_mut = false;
							Is_SC_APC_One_Hit_pt_mut = false;
							Is_SC_APC_Two_Hit_mut = false;
							Is_SC_BC_One_Hit_mut = false;
							Is_SC_BC_Two_Hit_mut = false;
							Is_BC_affected = false;

							for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
						 	cell_iter != cell_population.End();
						 	++cell_iter)
							{
								double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
								if (cell_height <= 0.5) 
								{
									cell_iter->SetMutationState(p_state_0); // all cells of the lowest TA cell row lose their mutations.
								}
							}
						}
						// end of copy
					}

					if (!NewStemCell && Is_SC_APC_One_Hit_pt_mut) // WT and only point mutation so far, second hit now
					{
						if (RandomNumberGenerator::Instance()->ranf() < 0.5) // LOH of WT allele
						{
							cout << "SC " << current_SC+1 << " undergoes APC-- due to LOH in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
							Is_SC_APC_One_Hit_pt_mut = false; 
							Is_SC_APC_Two_Hit_mut = true;
							ApcHitProb = 0.0;
						}
						else // LOH of mutated allele
						{	
							cout << "APC mutation is lost in SC " << current_SC+1 << " due to APC LOH in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event

							Is_SC_APC_One_Hit_pt_mut = false;
							Is_SC_APC_One_Hit_LOH_mut = true;
						}
					}

					if (!NewStemCell && !Is_SC_APC_One_Hit_LOH_mut && !Is_SC_APC_One_Hit_pt_mut && !Is_SC_APC_Two_Hit_mut) // first APC mutation
					{
						cout << "SC " << current_SC+1 << " undergoes APC+- due to LOH in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
						Is_SC_APC_One_Hit_LOH_mut = true;
						ApcHitProb *= 0.5;	
					}
				}

				if (Is_SC_APC_One_Hit_LOH_mut) // APC +- stem cell
				{
					for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				     cell_iter != cell_population.End();
				     ++cell_iter)
					{
				    	double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

				    	if (cell_height <= 0.5) // all cells of the lowest TA cell row are mutated. The new mutation status depends on previous mutations.
				    	{
							if (cell_iter->GetMutationState()->IsType<MMRTwoHitCellMutationState>()) 
							{
								cell_iter->SetMutationState(p_state_13); // APC+- and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<BCLOHMMRTwoHitCellMutationState>()) 
							{
								cell_iter->SetMutationState(p_state_20); // APC+-, CTNNB1 LOH and MMR
							}
				       		else if (cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_12); // APC+-
							}
							else if (cell_iter->GetMutationState()->IsType<BCLOHCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_19); // APC+- and CTNNB1 LOH
							}
							else if (cell_iter->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_14); // APC+- and beta-catenin
							}
							else if (cell_iter->GetMutationState()->IsType<MMRBCOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_15); // APC+-, beta-catenin and MMR
							}
				    	}
					}
				}

				// Now do point mutations.
				
				if (RandomNumberGenerator::Instance()->ranf() < ApcHitProb)
				{	
					if (Is_SC_APC_One_Hit_pt_mut || Is_SC_APC_One_Hit_LOH_mut ) // second APC hit
					{
						cout << "SC " << current_SC+1 << " undergoes APC-- due to point mutation in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
						Is_SC_APC_One_Hit_pt_mut = false;
						Is_SC_APC_Two_Hit_mut = true;
					}

					if (!Is_SC_APC_One_Hit_LOH_mut && !Is_SC_APC_One_Hit_pt_mut && !Is_SC_APC_Two_Hit_mut) // first APC mutation
					{
						cout << "SC " << current_SC+1 << " undergoes APC+- due to point mutation in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
						Is_SC_APC_One_Hit_pt_mut = true;
					}
				}

				if (Is_SC_APC_One_Hit_pt_mut) // APC +- stem cell
				{
					for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				     cell_iter != cell_population.End();
				     ++cell_iter)
					{
				    	double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

				    	if (cell_height <= 0.5) // all cells of the lowest TA cell row are mutated. The new mutation status depends on previous mutations.
				    	{
							if (cell_iter->GetMutationState()->IsType<MMRTwoHitCellMutationState>()) 
							{
								cell_iter->SetMutationState(p_state_5); // APC+- and MMR
							}
				       		else if (cell_iter->GetMutationState()->IsType<BCLOHMMRTwoHitCellMutationState>()) 
							{
								cell_iter->SetMutationState(p_state_18); // APC+-, CTNNB1 LOH and MMR
							}
				       		else if (cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_12); // APC+-
							}
							else if (cell_iter->GetMutationState()->IsType<BCLOHCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_17); // APC+- and CTNNB1 LOH
							}
							else if (cell_iter->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_8); // APC+- and beta-catenin
							}
							else if (cell_iter->GetMutationState()->IsType<MMRBCOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_9); // APC+-, beta-catenin and MMR
							}
				    	}
					}
				}

				if (Is_SC_APC_Two_Hit_mut) // APC -- stem cell. Whether due to LOH or point mutation does not make a difference.
				{
	
					for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				     cell_iter != cell_population.End();
				     ++cell_iter)
					{
				   		double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

				    	if (cell_height <= 0.5) // all cells of the lowest TA cell row are mutated. The new mutation status depends on previous mutations.
				    	{
						
							if (cell_iter->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<ApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>() || cell_iter->GetMutationState()->IsType<ApcLOHCellMutationState>()|| cell_iter->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>()) 
							{			
								cell_iter->SetMutationState(p_state_3); // APC-- (or beta-catenin and APC--, which is functionally equivalent in the cell cycle model)
							}
				       		
							else if (cell_iter->GetMutationState()->IsType<MMRTwoHitCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHMMRTwoHitCellMutationState>() || cell_iter->GetMutationState()->IsType<MMRApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHMMRApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<MMRApcLOHCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHMMRApcLOHCellMutationState>() || cell_iter->GetMutationState()->IsType<MMRBCOneHitApcOneHitCellMutationState>()  || cell_iter->GetMutationState()->IsType<MMRBCOneHitApcLOHCellMutationState>()|| cell_iter->GetMutationState()->IsType<MMRBCOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_6); // APC-- and MMR (or beta-catenin, APC-- and MMR, which is functionally equivalent)
							}
				    	}
					}
				}

				// Lastly, we do beta-catenin mutations. We start with LOH, which is only relevant for cell death, in case either another point mutation occurs, or if the cell becomes MLH1-deficient.

				if (RandomNumberGenerator::Instance()->ranf() < BetaCateninLOHProb)
				{	
					bool NewStemCell = false; // again track whether stem cell dies			
					
					if (Is_BC_affected) // second LOH event on the other allele (probability was halved)
					{
						NewStemCell = true; // cell dies

						// copy from above again:
						cout << "SC " << current_SC+1 << " of crypt " << j << " in cluster " << rank  <<  " is lost due to CTNNB1 LOH at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of SC loss
					
						// The stem cell loss is compensated for by symmetric division of an adjacent stem cell.				
							
						if (current_SC < SC_number)
						{
							all_mut_states[current_SC] = all_mut_states[current_SC+1];
							cout << "The loss is replenished by SC " << current_SC+2 << "." << endl;
						}
						else
						{
							all_mut_states[current_SC] = all_mut_states[0];
							cout << "The loss is replenished by SC " << 1 << "." << endl;
						}

						// A new stem cell is chosen for populating the crypt in the meantime.

						// Change cell numbering to make things easier, lost cell is now the "last" (i.e. the cell index) cell (in case it was not already).
						std::vector<bool> temp = all_mut_states[SC_index]; // store status of the current "last" cell in the ordering temporarily
						all_mut_states[SC_index] = all_mut_states[current_SC]; // first swap
						all_mut_states[current_SC]= temp; // second swap
						
						if (SC_index != current_SC)
						{
							cout << "We changed the cell numbering. SC " << current_SC+1 << " is now SC " << SC_index+1 << " and vice versa." << endl;
						}

						double p_old_SC = (double)SC_index/(double)SC_number; // the probability that the next stem cell is an "old" one

						if (RandomNumberGenerator::Instance()->ranf() < p_old_SC) // old stem cell
						{
							for (int i=1; i < SC_index; i++)
							{
								// replace mutation status by the ones of the old stem cell						
			
								if (RandomNumberGenerator::Instance()->ranf() < (double)(SC_index - i)/(double)(SC_index +1 -i))
								{
									Is_SC_MMR_mut = all_mut_states[i][0];
									Is_SC_MMR_LOH_mut = all_mut_states[i][1];
									Is_SC_MLH1_LOH_mut = all_mut_states[i][2];
									Is_SC_APC_One_Hit_LOH_mut = all_mut_states[i][3];
									Is_SC_APC_One_Hit_pt_mut = all_mut_states[i][4];
									Is_SC_APC_Two_Hit_mut = all_mut_states[i][5];
									Is_SC_BC_One_Hit_mut = all_mut_states[i][6];
									Is_SC_BC_Two_Hit_mut = all_mut_states[i][7];
									Is_BC_affected = all_mut_states[i][8];
									if (i == (SC_index-1)) // only last one is relevant, i.e. the "previous" cell
									{
										cout << "It is now populated by SC " << i+1 << " again." << endl;
										current_SC = i;										
									}
								}
								else
								{
									Is_SC_MMR_mut = all_mut_states[i-1][0];
									Is_SC_MMR_LOH_mut = all_mut_states[i-1][1];
									Is_SC_MLH1_LOH_mut = all_mut_states[i-1][2];
									Is_SC_APC_One_Hit_LOH_mut = all_mut_states[i-1][3];
									Is_SC_APC_One_Hit_pt_mut = all_mut_states[i-1][4];
									Is_SC_APC_Two_Hit_mut = all_mut_states[i-1][5];
									Is_SC_BC_One_Hit_mut = all_mut_states[i-1][6];
									Is_SC_BC_Two_Hit_mut = all_mut_states[i-1][7];
									Is_BC_affected = all_mut_states[i-1][8];
									cout << "It is now populated by SC " << i << " again." << endl;
									current_SC = i-1;
									break;
								}
							}
						}
						else // new stem cell
						{
							SC_index++;
							current_SC = SC_index;
							cout << "It is now populated by a new SC, namely SC " << current_SC+1 << ". Mutations are lost." << endl; 
						
							//Reset boolean variables
							Is_SC_MMR_mut = false;
							Is_SC_MMR_LOH_mut = false;
							Is_SC_MLH1_LOH_mut = false;
							Is_SC_APC_One_Hit_LOH_mut = false;
							Is_SC_APC_One_Hit_pt_mut = false;
							Is_SC_APC_Two_Hit_mut = false;
							Is_SC_BC_One_Hit_mut = false;
							Is_SC_BC_Two_Hit_mut = false;
							Is_BC_affected = false;

							for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
						 	cell_iter != cell_population.End();
						 	++cell_iter)
							{
								double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
								if (cell_height <= 0.5) 
								{
									cell_iter->SetMutationState(p_state_0); // all cells of the lowest TA cell row lose their mutations.
								}
							}
						}
						// end of copy
					}
					
					if (!NewStemCell && Is_SC_BC_One_Hit_mut) // WT and only point mutation so far, second hit now
					{
						if (RandomNumberGenerator::Instance()->ranf() < 0.5) // LOH of WT allele
						{
							cout << "SC " << current_SC+1 << " undergoes CTNNB1++ due to LOH in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
							Is_SC_BC_One_Hit_mut = false;
							Is_BC_affected = true;
							Is_SC_BC_Two_Hit_mut = true;
							BetaCateninProb = 0.0;
						}
						else // LOH of mutated allele
						{	
							cout << "CTNNB1 mutation is lost in SC " << current_SC+1 << " due to CTNNB1 LOH in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event

							Is_SC_BC_One_Hit_mut = false;
							Is_BC_affected = true;
						}
					}

					if (!NewStemCell && !Is_BC_affected && !Is_SC_BC_One_Hit_mut && !Is_SC_BC_Two_Hit_mut) // first CTNNB1 mutation
					{
						cout << "SC " << current_SC+1 << " undergoes CTNNB1 LOH in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
						Is_BC_affected = true;
						BetaCateninProb *= 0.5;	
					}
					
					// Now take care of the MLH1 LOH

					if (!NewStemCell && MLH1 && !Is_SC_MLH1_LOH_mut && RandomNumberGenerator::Instance()->ranf() < 0.8) // concurrent MLH1 LOH
					{	

						cout << "SC " << current_SC+1 << " undergoes MLH1 LOH due to CTNNB1 LOH in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output time of event

						if (Is_SC_MMR_LOH_mut && RandomNumberGenerator::Instance()->ranf() < 0.5) // one LOH so far, now loss of WT/pt mut allele
						{
							Is_SC_MLH1_LOH_mut = true;
							Is_SC_MMR_mut = true; // cell is now MMR-deficient (or it was already)
						}

						if (!Is_SC_MMR_LOH_mut && Is_SC_MMR_mut) // no LOH but MMR-deficient, i.e. two pt muts
						{
							Is_SC_MMR_LOH_mut = true;
						}

						if (!Is_SC_MMR_LOH_mut && !Is_SC_MMR_mut) // only one pt mut ("wild-type")
						{
							if (RandomNumberGenerator::Instance()->ranf() < 0.5) // LOH of mutated allele
							{
								Is_SC_MMR_LOH_mut = true;
							}
							else // LOH of WT allele
							{
								Is_SC_MMR_LOH_mut = true;
								Is_SC_MMR_mut = true;
							}
						}
			
						if (Is_SC_MMR_mut) // cell is now MMR-deficient (or already was), mutate lowest cell row as above
						{	
				
							for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
						 	cell_iter != cell_population.End();
						 	++cell_iter)
							{
								double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

								if (cell_height <= 0.5) // all cells of the lowest TA cell row are mutated. The new mutation status depends on previous mutations.
								{
									if (cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_1); // MMR mutation
									}
									else if (cell_iter->GetMutationState()->IsType<BCLOHCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_16); // CTNNB1 LOH and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<ApcOneHitCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_5); // APC+- and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<ApcLOHCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_13); 
									}
									else if (cell_iter->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_18); // APC+-, CTNNB1 LOH and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_20);
									}
									else if (cell_iter->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_6); // APC-- and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_7); // beta-catenin + and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<BetaCateninTwoHitCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_11); // beta-catenin ++ and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_9); // beta-catenin +, APC+- and MMR
									}
									else if (cell_iter->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
									{
										cell_iter->SetMutationState(p_state_15); 
									}
								}
							}
						}
					}
				}

				// Now do point mutations.
	
				if (RandomNumberGenerator::Instance()->ranf() < BetaCateninProb)
				{	

					if (Is_SC_BC_One_Hit_mut || Is_BC_affected) // second beta-catenin mutation
					{
						cout << "SC " << current_SC+1 << " undergoes beta-catenin ++ mutation in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
						Is_SC_BC_One_Hit_mut = false;
						Is_SC_BC_Two_Hit_mut = true; 
					}

					if (!Is_SC_BC_One_Hit_mut && !Is_SC_BC_Two_Hit_mut && !Is_BC_affected) // first beta-catenin mutation
					{
						
						cout << "SC " << current_SC+1 << " undergoes beta-catenin + mutation in crypt " << j << " in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
						Is_SC_BC_One_Hit_mut = true;
					}
				}

				if (Is_SC_BC_One_Hit_mut)
				{
					for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
		        	 cell_iter != cell_population.End();
		        	 ++cell_iter)
		    		{
		        		double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

		        		if (cell_height <= 0.5) // all cells of the lowest TA cell row are mutated. The new mutation status depends on previous mutations.
		        		{
							if (cell_iter->GetMutationState()->IsType<MMRTwoHitCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHMMRTwoHitCellMutationState>()) 
							{
								cell_iter->SetMutationState(p_state_7); // beta-catenin and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<ApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>()) 
							{			
								cell_iter->SetMutationState(p_state_8); // beta-catenin and APC+-
							}
							else if (cell_iter->GetMutationState()->IsType<ApcLOHCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>()) 
							{			
								cell_iter->SetMutationState(p_state_14); 
							}
		           			else if (cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHCellMutationState>() )
							{
								cell_iter->SetMutationState(p_state_4); // beta-catenin mutation
							}
							else if (cell_iter->GetMutationState()->IsType<MMRApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHMMRApcOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_9); // beta-catenin, APC+- and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<MMRApcLOHCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHMMRApcLOHCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_15); 
							}
		        		}
					}
				}

				if (Is_BC_affected)
				{
					for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
		        	 cell_iter != cell_population.End();
		        	 ++cell_iter)
		    		{
		        		double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

		        		if (cell_height <= 0.5) // all cells of the lowest TA cell row are mutated. The new mutation status depends on previous mutations.
		        		{
							if (cell_iter->GetMutationState()->IsType<MMRTwoHitCellMutationState>() || cell_iter->GetMutationState()->IsType<MMRBCOneHitCellMutationState>()) 
							{
								cell_iter->SetMutationState(p_state_16); // CTNNB1 LOH and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<ApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>()) 
							{			
								cell_iter->SetMutationState(p_state_17); // CTNNB1 LOH and APC+-
							}
							else if (cell_iter->GetMutationState()->IsType<ApcLOHCellMutationState>() || cell_iter->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>()) 
							{			
								cell_iter->SetMutationState(p_state_19); 
							}
		           			else if (cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>() || cell_iter->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>() )
							{
								cell_iter->SetMutationState(p_state_21); // CTNNB1 LOH
							}
							else if (cell_iter->GetMutationState()->IsType<MMRApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<MMRBCOneHitApcOneHitCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_18); // CTNNB1 LOH APC+- and MMR
							}
							else if (cell_iter->GetMutationState()->IsType<MMRApcLOHCellMutationState>() || cell_iter->GetMutationState()->IsType<MMRBCOneHitApcLOHCellMutationState>())
							{
								cell_iter->SetMutationState(p_state_20); 
							}
		        		}
					}
				}


				if (Is_SC_BC_Two_Hit_mut)
				{

					for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
		        	 cell_iter != cell_population.End();
		        	 ++cell_iter)
		    		{
		        		double cell_height = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

		        		if (cell_height <= 0.5) // all cells of the lowest TA cell row are mutated. The new mutation status depends on previous mutations.
		        		{
							if (cell_iter->GetMutationState()->IsType<MMRTwoHitCellMutationState>() || cell_iter->GetMutationState()->IsType<MMRApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<MMRApcLOHCellMutationState>() || cell_iter->GetMutationState()->IsType<MMRBCOneHitCellMutationState>() ||  cell_iter->GetMutationState()->IsType<MMRBCOneHitApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<MMRBCOneHitApcLOHCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHMMRTwoHitCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHMMRApcOneHitCellMutationState>() || cell_iter->GetMutationState()->IsType<BCLOHMMRApcLOHCellMutationState>()) // all MMR mutations
							{
								cell_iter->SetMutationState(p_state_11); // beta-catenin++ and MMR
							}
							else // the rest
							{			
								cell_iter->SetMutationState(p_state_10); // beta-catenin ++
							}
		        		}
					}
				}
			}

			// This marks the end of the crypt simulation. We now check whether the crypt is MMR-deficient and output how many crypt fissions have occurred.
		
			double mmr_percentage = (double)cell_population.GetNumMMRCells()/(double)cell_population.GetNumRealCells(); // once again check for MMR-deficiency at the end
			if (mmr_percentage > 0.8)
			{
				cout << "Crypt " << j  << " in cluster " << rank << " gives rise to a MMR-DCF consisting of " << p_killer2->GetCryptCounter() << " crypts." << endl;
				MMR_crypt_no += p_killer2->GetCryptCounter(); // increase number of MMR-def. crypts
				MMR_DCF_no++;  // increment number of crypt foci
			}
		
			total_crypt_no += p_killer2->GetCryptCounter(); // increase total crypt number

		   /* Finally, we must tidy up by destroying the {{{WntConcentration}}} and {{{SimulationTime}}}
		    * singleton objects. This avoids memory leaks occurring.
			* We must also reset the start time to 0 to be able to start the next simulation. */

		   WntConcentration<2>::Destroy();
		   SimulationTime::Destroy();
		   SimulationTime::Instance()->SetStartTime(0.0);
		   RandomNumberGenerator::Instance()->Reseed(getpid()); // reseed for the next crypt 
		
		}
	
	// This marks the end of the outer for-loop, i.e. the cluster simulation. We first calculate the number of new crypts..

	int new_crypts = total_crypt_no - initial_crypt_no; 
	
	// .. and summarize the results of the simulation.
	
	cout << "End of simulation in cluster " << rank << "." << endl; 
	cout << "The total number of crypts in cluster " << rank << " is " << total_crypt_no << ", which is an increase of " << new_crypts << "." << endl;
	cout << "The total number of MMR-deficient crypts in cluster " << rank << " is " << MMR_crypt_no << "." << endl;
	cout << "The number of MMR-DCF in cluster " << rank << " is " << MMR_DCF_no << "." << endl;

	MPI_Finalize();

	}

};

#endif /*TESTCOLONICCRYPTSIMULATION_HPP_*/
