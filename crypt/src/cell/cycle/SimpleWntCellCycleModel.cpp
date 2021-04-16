/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
/* The first two header files are needed for parallelization.*/
#include <mpi.h>
#include <stdio.h>

#include "SimpleWntCellCycleModel.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "BetaCateninTwoHitCellMutationState.hpp"
#include "MMRTwoHitCellMutationState.hpp"
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

SimpleWntCellCycleModel::SimpleWntCellCycleModel()
    : mUseCellProliferativeTypeDependentG1Duration(false),
      mWntStemThreshold(0.9), // do not need this
      mWntTransitThreshold(0.75),
      mWntLabelledThreshold(0.75)
{
}

SimpleWntCellCycleModel::SimpleWntCellCycleModel(const SimpleWntCellCycleModel& rModel)
   : AbstractSimplePhaseBasedCellCycleModel(rModel),
     mUseCellProliferativeTypeDependentG1Duration(rModel.mUseCellProliferativeTypeDependentG1Duration),
     mWntStemThreshold(rModel.mWntStemThreshold),
     mWntTransitThreshold(rModel.mWntTransitThreshold),
     mWntLabelledThreshold(rModel.mWntLabelledThreshold)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration and the cell's proliferative type are
     * (re)set as soon as InitialiseDaughterCell() is called on the
     * new cell-cycle model.
     */
}

AbstractCellCycleModel* SimpleWntCellCycleModel::CreateCellCycleModel()
{
    return new SimpleWntCellCycleModel(*this);
}

bool SimpleWntCellCycleModel::GetUseCellProliferativeTypeDependentG1Duration() const
{
    return mUseCellProliferativeTypeDependentG1Duration;
}

void SimpleWntCellCycleModel::SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration)
{
    mUseCellProliferativeTypeDependentG1Duration = useCellProliferativeTypeDependentG1Duration;
}

void SimpleWntCellCycleModel::SetG1Duration()
{
    assert(mpCell != nullptr);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
	// do not explicitly consider stem cells anymore
	if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
	{
	    if (mUseCellProliferativeTypeDependentG1Duration)
	    {
	        mG1Duration = p_gen->NormalRandomDeviate(GetStemCellG1Duration(), 0.5);
	    }
	    else
	    {
	        // Normally stem cells should behave just like transit cells in a Wnt simulation
	        mG1Duration = p_gen->NormalRandomDeviate(GetTransitCellG1Duration(), 0.5);
	    }
	}
	else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
	{
	    mG1Duration = p_gen->NormalRandomDeviate(GetTransitCellG1Duration(), 0.5);
	}
	else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
	{
	    mG1Duration = DBL_MAX;
	}
	else
	{
	    NEVER_REACHED;
	}

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < mMinimumGapDuration)
    {
        mG1Duration = mMinimumGapDuration;
    }
}

double SimpleWntCellCycleModel::GetWntLevel() const
{
    assert(mpCell != nullptr);
    double level = 0;

    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        default:
            NEVER_REACHED;
    }
    return level;
}

WntConcentrationType SimpleWntCellCycleModel::GetWntType()
{
    WntConcentrationType wnt_type;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case UNSIGNED_UNSET:
        {
            // If you trip this you have tried to use a simulation without setting the dimension.
            NEVER_REACHED;
        }
        default:
            NEVER_REACHED;
    }
    return wnt_type;
}

void SimpleWntCellCycleModel::UpdateCellCyclePhase()
{

    // The cell can divide if the Wnt concentration >= wnt_division_threshold

    double wnt_division_threshold = DBL_MAX;

    // Set up under what level of Wnt stimulus a cell will divide

    if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>() || mpCell->GetMutationState()->IsType<MMRTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHCellMutationState>())
    {
        wnt_division_threshold = mWntTransitThreshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<ApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<MMRApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRApcLOHCellMutationState>())
    {
        // should be less than healthy values, all APC+- dominated mutations
        wnt_division_threshold = 0.8*mWntTransitThreshold;
    }
    else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitCellMutationState>())
    {
        // less than above value, all beta-catenin-dominated mutations
        wnt_division_threshold = 0.9*mWntTransitThreshold;
	}
	else if (mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcLOHCellMutationState>())
	{
		wnt_division_threshold = 0.72*mWntTransitThreshold; // product of the above two factors
	}

    else if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRApcTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<BetaCateninTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCTwoHitCellMutationState>())
    {
        // set to zero (no Wnt-dependence) for all APC-- dominated mutations and beta-catenin double hits
        wnt_division_threshold = 0.0;
    }
    else
    {
        NEVER_REACHED;
    }

    if (mpCell->HasCellProperty<CellLabel>())
    {
        wnt_division_threshold = mWntLabelledThreshold;
    }

    double wnt_level = GetWntLevel();


    // Set the cell type to TransitCellProliferativeType if the Wnt stimulus exceeds wnt_division_threshold. We do not explicitly consider stem cells anymore.
	/*
     * This method is usually called within a CellBasedSimulation, after the CellPopulation
     * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
     * CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
     * would be creating a new CellPropertyRegistry. In this case the cell proliferative
     * type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
     * would be incorrect. We must therefore access the CellProliferativeType via the cell's
     * CellPropertyCollection.
     */

	// Cell is a TA cell if Wnt level is above the particular threshold.
    if (wnt_level >= wnt_division_threshold)
    {
   	    boost::shared_ptr<AbstractCellProperty> p_transit_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);
    }
    else
    {
        // The cell is set to be an FD cell and hence is in G0 phase.
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
    AbstractSimplePhaseBasedCellCycleModel::UpdateCellCyclePhase();
}


	// We create shared pointers to the mutation states we wish to use.

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



void SimpleWntCellCycleModel::InitialiseDaughterCell()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // initialize rank to put out in which cluster the mutation occurs.

	// We define several boolean variables which facilitate changing the mutation rates and the ouput of mutations.

	bool Is_Cell_APC_One_Hit_pt_mut = false;
	bool Is_Cell_APC_One_Hit_LOH_mut = false;
	bool Is_Cell_APC_Two_Hit_mut = false;
	bool Is_Cell_MMR_mut = false;
	bool Is_Cell_BC_One_Hit_mut = false;
	bool Is_Cell_BC_Two_Hit_mut = false;
	bool Is_BC_affected = false; // LOH of CTNNB1

	// We first implement mutation-induced cell death.

	double mutInducedDeathRate = 0.0001; // death rate caused by mutations, corresponds to about 1/30 % deleterious mutations


	if (mpCell->GetMutationState()->IsType<MMRTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<MMRApcTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRApcLOHCellMutationState>())
	{
		mutInducedDeathRate *= 100.0; // Increase mutation-induced death probability of MMR-mutant cells hundred-fold.
		Is_Cell_MMR_mut = true;
	}


	if (RandomNumberGenerator::Instance()->ranf() < mutInducedDeathRate)
	{
		mpCell->StartApoptosis(); // cell dies shortly after division
	}


	/* We need to distinguish between MLH1 and MSH2 germline mutation carriers, since MLH1 LOH of the WT allele can induce a CTNNB1 double hit, or lead to the loss of a CTNNB1 mutation.
	 * To simulate different patients, we need to alter the boolean variable here and in TestColonicCryptSimulation.hpp. */

	bool MLH1 = false;

	// Set up the random mutations..

	double MMRMutProb; // Probability of a sporadic MMR mutation.

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
	double BetaCateninProb = 0.0000000375; // Probability of a sporadic beta-catenin mutation.
	double BetaCateninLOHProb = 0.000010116; // Probability of Beta-Catenin LOH (important if second allele becomes mutated)

	// We increase the mutation probability of point mutations (as we do with death) 100-fold for MMR-deficient cells.

	if (Is_Cell_MMR_mut)
	{
		ApcHitProb *= 100.0;
		BetaCateninProb *= 100.0;
	}


	if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRApcOneHitCellMutationState>())
	{
		ApcHitProb *= 0.5;
		Is_Cell_APC_One_Hit_pt_mut = true;
	}

	// If the first hit was due to LOH, the probability of another LOH event is halved.

	if (mpCell->GetMutationState()->IsType<ApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<MMRApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRApcLOHCellMutationState>())
	{
		ApcHitProb *= 0.5;
		ApcLOHProb *= 0.5;
		Is_Cell_APC_One_Hit_LOH_mut = true;
	}

	// mark APC -- cells

	if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRApcTwoHitCellMutationState>())
	{
		Is_Cell_APC_Two_Hit_mut = true;
	}

	// mark beta-catenin single hits and halve the probability

	if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>()  || mpCell->GetMutationState()->IsType<MMRBCOneHitApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>() ||  mpCell->GetMutationState()->IsType<MMRBCOneHitApcLOHCellMutationState>())
	{
		BetaCateninProb *= 0.5;
		Is_Cell_BC_One_Hit_mut = true;
	}

	// mark beta-catenin LOHs

	if (mpCell->GetMutationState()->IsType<BCLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>()  || mpCell->GetMutationState()->IsType<BCLOHMMRApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>() ||  mpCell->GetMutationState()->IsType<BCLOHMMRApcLOHCellMutationState>())
	{
		Is_BC_affected = true;
		BetaCateninProb *= 0.5;
		BetaCateninLOHProb *= 0.5;
	}

	// mark beta-catenin double hits

	if (mpCell->GetMutationState()->IsType<BetaCateninTwoHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCTwoHitCellMutationState>())
	{
		Is_Cell_BC_Two_Hit_mut = true;
	}



	// We include random mutations upon cell division here. We broadly perform the same procedure as the SC mutation in the main code TestColonicCryptSimulation.hpp.

	// MMR point mutation:
	if (RandomNumberGenerator::Instance()->ranf() < MMRMutProb)
	{

		if (!Is_Cell_MMR_mut) // novel MMR mutation
		{
			cout << "TA cell MMR point mutation in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
		}

		if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
		{
			mpCell->SetMutationState(p_state_1); // MMR mutation
		}
		else if (mpCell->GetMutationState()->IsType<BCLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_16);
		}
		else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_5); // APC+- and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_18);
		}
		else if (mpCell->GetMutationState()->IsType<ApcLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_13);
		}
		else if (mpCell->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_20);
		}
		else if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_6); // APC-- and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_7); // beta-catenin + and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BetaCateninTwoHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_11); // beta-catenin ++ and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_9); // beta-catenin +, APC+- and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_15);
		}
	}

	// MMR LOH:

	if (RandomNumberGenerator::Instance()->ranf() < MMRLOHProb)
	{

		// Output the time of mutation in all cases. Even if the allele already has a point mutation, the LOH might still be relevant.

		cout << "TA cell MMR LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl;

		// First, take care of concurrent CTNNB1 LOH

		if (MLH1 && RandomNumberGenerator::Instance()->ranf() < 0.8) // MLH1 LOH affects CTNNB1 in 80% of cases.
		{
			if (!Is_BC_affected) // no CTNNB1 LOH so far
			{
				Is_BC_affected = true; // BC is lost.
				BetaCateninProb *= 0.5; // There is one less WT allele to be mutated, unless a CTNNB1 mutation was present on the lost allele (see below).
				BetaCateninLOHProb *= 0.5;

				// Change the mutation state accordingly

				if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
				{
					mpCell->SetMutationState(p_state_21); // CTNNB1 LOH
				}
				else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
				{
					mpCell->SetMutationState(p_state_17); // APC+- and CTNNB1 LOH
				}
				else if (mpCell->GetMutationState()->IsType<ApcLOHCellMutationState>())
				{
					mpCell->SetMutationState(p_state_19);
				}
				else if (mpCell->GetMutationState()->IsType<MMRTwoHitCellMutationState>())
				{
					mpCell->SetMutationState(p_state_16); // CTNNB1 LOH and MMR
				}
				else if (mpCell->GetMutationState()->IsType<MMRApcOneHitCellMutationState>())
				{
					mpCell->SetMutationState(p_state_18); // CTNNB1 LOH, APC+- and MMR
				}
				else if (mpCell->GetMutationState()->IsType<MMRApcLOHCellMutationState>())
				{
					mpCell->SetMutationState(p_state_20);
				}
				else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
				{
					if (RandomNumberGenerator::Instance()->ranf() < 0.5) // BC mutation is on other allele
					{
						mpCell->SetMutationState(p_state_10); // beta-catenin ++
						cout << "TA cell beta-catenin double hit due to MLH1 LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
					}
					else // BC mutation is on lost allele
					{
						BetaCateninProb *= 2.0; // correct factor from above, probability was already halved
						mpCell->SetMutationState(p_state_21); // CTNNB1 mut is lost
						cout << "TA cell beta-catenin mutation is lost due to MLH1 LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
					}
				}
				else if (mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>())
				{

					if (RandomNumberGenerator::Instance()->ranf() < 0.5) // BC mutation is on other allele
					{
						mpCell->SetMutationState(p_state_10); // beta-catenin ++ (and APC+-)
						cout << "TA cell beta-catenin double hit due to MLH1 LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
					}
					else // BC mutation is on lost allele
					{
						BetaCateninProb *= 2.0; // correct factor from above, probability was already halved
						mpCell->SetMutationState(p_state_17); // APC +-, beta-catenin mut is lost
						cout << "TA cell beta-catenin mutation is lost due to MMR LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
					}
				}
				else if (mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
				{
					if (RandomNumberGenerator::Instance()->ranf() < 0.5) // BC mutation is on other allele
					{
						mpCell->SetMutationState(p_state_10); // beta-catenin ++ (and APC+-)
						cout << "TA cell beta-catenin double hit due to MMR LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
					}
					else // BC mutation is on lost allele
					{
						BetaCateninProb *= 2.0; // correct factor from above, probability was already halved
						mpCell->SetMutationState(p_state_19); // APC +-, beta-catenin mut is lost
						cout << "TA cell beta-catenin mutation is lost due to MMR LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
					}
				}
			}
			else
			{
				if (RandomNumberGenerator::Instance()->ranf() < 0.5)
				{
					mpCell->StartApoptosis(); // cell dies due to two LOH events. Otherwise, nothing changes.
				}
			}
		}

		// Now take care of the MMR mutation

		if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
		{
			mpCell->SetMutationState(p_state_1); // MMR mutation
		}
		else if (mpCell->GetMutationState()->IsType<BCLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_16);
		}
		else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_5); // APC+- and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_18);
		}
		else if (mpCell->GetMutationState()->IsType<ApcLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_13);
		}
		else if (mpCell->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_20);
		}
		else if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_6); // APC-- and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_7); // beta-catenin + and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BetaCateninTwoHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_11); // beta-catenin ++ and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_9); // beta-catenin +, APC +- and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_15); // beta-catenin +, APC +- and MMR
		}
	}

	// APC point mutation:

	if (RandomNumberGenerator::Instance()->ranf() < ApcHitProb)
	{

		if (Is_Cell_APC_One_Hit_pt_mut || Is_Cell_APC_One_Hit_LOH_mut) // second hit
		{
			cout << "TA cell APC-- point mutation in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
		}

		if (!Is_Cell_APC_One_Hit_pt_mut && !Is_Cell_APC_One_Hit_LOH_mut && !Is_Cell_APC_Two_Hit_mut) // novel APC mutation.
		{
			cout << "TA cell APC+- point mutation in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
			Is_Cell_APC_One_Hit_pt_mut = true;
			ApcLOHProb *= 0.5;
		}

		if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
		{
			mpCell->SetMutationState(p_state_2); // APC+-
		}
		else if (mpCell->GetMutationState()->IsType<BCLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_17); // APC+- and CTNNB1 LOH
		}
		else if (mpCell->GetMutationState()->IsType<MMRTwoHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_5); // APC+- and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BCLOHMMRTwoHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_18); // APC+-, CTNNB1 LOH and MMR
		}
		else if (mpCell->GetMutationState()->IsType<MMRBCOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_9); // APC+-, beta-catenin+ and MMR
		}
		else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<ApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_3); // APC-- (or beta-catenin+ and APC--, which is functionally equivalent in the cell cycle model)
		}
		else if (mpCell->GetMutationState()->IsType<MMRApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_6); // APC-- and MMR (or beta-catenin+, APC-- and MMR, which is functionally equivalent)
		}
		else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_8); // APC+- and beta-catenin+
		}
	}

	// APC LOH:

	if (RandomNumberGenerator::Instance()->ranf() < ApcLOHProb)
	{

		if (!Is_Cell_APC_One_Hit_pt_mut && !Is_Cell_APC_One_Hit_LOH_mut && !Is_Cell_APC_Two_Hit_mut) // novel LOH
		{
			cout << "TA cell APC+- LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
		}

		if (Is_Cell_APC_One_Hit_pt_mut) // second hit
		{
			cout << "TA cell APC-- LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
		}

		if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
		{
			mpCell->SetMutationState(p_state_12); // APC+-
		}
		else if (mpCell->GetMutationState()->IsType<BCLOHCellMutationState>())
		{
			mpCell->SetMutationState(p_state_19); // APC+- and CTNNB1 LOH
		}
		else if (mpCell->GetMutationState()->IsType<MMRTwoHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_13); // APC+- and MMR
		}
		else if (mpCell->GetMutationState()->IsType<BCLOHMMRTwoHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_20); // APC+-, CTNNB1 LOH and MMR
		}
		else if (mpCell->GetMutationState()->IsType<MMRBCOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_15); // APC+-, beta-catenin+ and MMR
		}
		else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>() ||  mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_3); // APC-- (or beta-catenin+ and APC--, which is functionally equivalent in the cell cycle model)
		}
		else if (mpCell->GetMutationState()->IsType<ApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
		{
			mpCell->StartApoptosis(); // two APC LOHs
		}
		else if (mpCell->GetMutationState()->IsType<MMRApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRApcOneHitCellMutationState>() ||  mpCell->GetMutationState()->IsType<MMRBCOneHitApcOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_6); // APC-- and MMR (or beta-catenin+, APC-- and MMR, which is functionally equivalent)
		}
		else if (mpCell->GetMutationState()->IsType<MMRApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<BCLOHMMRApcLOHCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcLOHCellMutationState>())
		{
			mpCell->StartApoptosis(); // two APC LOHs
		}
		else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
		{
			mpCell->SetMutationState(p_state_14); // APC+- and beta-catenin+
		}
	}

	// Beta-Catenin point mutation:
	if (RandomNumberGenerator::Instance()->ranf() < BetaCateninProb)
	{
		if (!Is_Cell_BC_One_Hit_mut && !Is_Cell_BC_Two_Hit_mut) // novel BC mutation
		{
			cout << "TA cell beta-catenin+ mutation in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
		}

		if (Is_Cell_BC_One_Hit_mut) // second hit
		{
			cout << "TA cell beta-catenin++ mutation in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output the time of mutation
		}

		if (Is_BC_affected) // only one beta-catenin allele left
		{
			if (Is_Cell_MMR_mut)
			{
				mpCell->SetMutationState(p_state_11); // beta-catenin++ and MMR
			}
			else
			{
				mpCell->SetMutationState(p_state_10); // beta-catenin++
			}
		}
		else	// normal dynamics
		{
			if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
			{
				mpCell->SetMutationState(p_state_4); // beta-catenin+ mutation
			}
			else if (mpCell->GetMutationState()->IsType<MMRTwoHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_7); // beta-catenin+ and MMR
			}
			else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_8); // beta-catenin+ and APC+-
			}
			else if (mpCell->GetMutationState()->IsType<ApcLOHCellMutationState>())
			{
				mpCell->SetMutationState(p_state_14);
			}
			else if (mpCell->GetMutationState()->IsType<MMRApcOneHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_9); // beta-catenin+, APC+- and MMR
			}
			else if (mpCell->GetMutationState()->IsType<MMRApcLOHCellMutationState>())
			{
				mpCell->SetMutationState(p_state_15);
			}
			else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
			{
				mpCell->SetMutationState(p_state_10); // beta-catenin++
			}
			else if (mpCell->GetMutationState()->IsType<MMRBCOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcOneHitCellMutationState>() || mpCell->GetMutationState()->IsType<MMRBCOneHitApcLOHCellMutationState>())
			{
				mpCell->SetMutationState(p_state_11); // MMR and beta-catenin++ (APC mutation has no influence anymore)
			}
		}
	}

	// Beta-Catenin LOH:

	if (RandomNumberGenerator::Instance()->ranf() < BetaCateninLOHProb)
	{
		if (Is_BC_affected)
		{
			mpCell->StartApoptosis(); // two LOHs
		}
		else
		{
			if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
			{
				mpCell->SetMutationState(p_state_21); // CTNNB1 LOH
			}
			else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_17); // APC+- and CTNNB1 LOH
			}
			else if (mpCell->GetMutationState()->IsType<ApcLOHCellMutationState>())
			{
			mpCell->SetMutationState(p_state_19);
			}
			else if (mpCell->GetMutationState()->IsType<MMRTwoHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_16); // CTNNB1 LOH and MMR
			}
			else if (mpCell->GetMutationState()->IsType<MMRApcOneHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_18); // CTNNB1 LOH, APC+- and MMR
			}
			else if (mpCell->GetMutationState()->IsType<MMRApcLOHCellMutationState>())
			{
				mpCell->SetMutationState(p_state_20);
			}
			else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
			{
				if (RandomNumberGenerator::Instance()->ranf() < 0.5) // BC mutation is on other allele
				{
					mpCell->SetMutationState(p_state_10); // beta-catenin ++
					cout << "TA cell beta-catenin double hit due to MLH1 LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
				}
				else // BC mutation is on lost allele
				{
					BetaCateninProb *= 2.0; // correct factor from above, probability was already halved
					mpCell->SetMutationState(p_state_21); // CTNNB1 mut is lost
					cout << "TA cell beta-catenin mutation is lost due to MLH1 LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
				}
			}
			else if (mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>())
			{

				if (RandomNumberGenerator::Instance()->ranf() < 0.5) // BC mutation is on other allele
				{
					mpCell->SetMutationState(p_state_10); // beta-catenin ++ (and APC+-)
					cout << "TA cell beta-catenin double hit due to MLH1 LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
				}
				else // BC mutation is on lost allele
				{
					BetaCateninProb *= 2.0; // correct factor from above, probability was already halved
					mpCell->SetMutationState(p_state_17); // APC +-, beta-catenin mut is lost
					cout << "TA cell beta-catenin mutation is lost due to MLH1 LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
				}
			}
			else if (mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
			{
				if (RandomNumberGenerator::Instance()->ranf() < 0.5) // BC mutation is on other allele
				{
					mpCell->SetMutationState(p_state_10); // beta-catenin ++ (and APC+-)
					cout << "TA cell beta-catenin double hit due to MLH1 LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
				}
				else // BC mutation is on lost allele
				{
					BetaCateninProb *= 2.0; // correct factor from above, probability was already halved
					mpCell->SetMutationState(p_state_19); // APC +-, beta-catenin mut is lost
					cout << "TA cell beta-catenin mutation is lost due to MLH1 LOH in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl; // output event
				}
			}
		}

		if (MLH1 && RandomNumberGenerator::Instance()->ranf() < 0.8) // CTNNB1 LOH affects MLH1 in 80% of cases.
		{
			if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
			{
				mpCell->SetMutationState(p_state_1); // MMR mutation
			}
			else if (mpCell->GetMutationState()->IsType<BCLOHCellMutationState>())
			{
				mpCell->SetMutationState(p_state_16);
			}
			else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_5); // APC+- and MMR
			}
			else if (mpCell->GetMutationState()->IsType<BCLOHApcOneHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_18);
			}
			else if (mpCell->GetMutationState()->IsType<ApcLOHCellMutationState>())
			{
				mpCell->SetMutationState(p_state_13);
			}
			else if (mpCell->GetMutationState()->IsType<BCLOHApcLOHCellMutationState>())
			{
				mpCell->SetMutationState(p_state_20);
			}
			else if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_6); // APC-- and MMR
			}
			else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_7); // beta-catenin + and MMR
			}
			else if (mpCell->GetMutationState()->IsType<BetaCateninTwoHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_11); // beta-catenin ++ and MMR
			}
			else if (mpCell->GetMutationState()->IsType<BCOneHitApcOneHitCellMutationState>())
			{
				mpCell->SetMutationState(p_state_9); // beta-catenin +, APC +- and MMR
			}
			else if (mpCell->GetMutationState()->IsType<BCOneHitApcLOHCellMutationState>())
			{
				mpCell->SetMutationState(p_state_15); // beta-catenin +, APC +- and MMR
			}

		}
	}

    AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();
}





bool SimpleWntCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

double SimpleWntCellCycleModel::GetWntStemThreshold() const
{
    return mWntStemThreshold;
}

void SimpleWntCellCycleModel::SetWntStemThreshold(double wntStemThreshold)
{
    assert(wntStemThreshold <= 1.0);
    assert(wntStemThreshold >= 0.0);
    mWntStemThreshold = wntStemThreshold;
}

double SimpleWntCellCycleModel::GetWntTransitThreshold() const
{
    return mWntTransitThreshold;
}

void SimpleWntCellCycleModel::SetWntTransitThreshold(double wntTransitThreshold)
{
    //assert(wntTransitThreshold <= 1.0);
    //assert(wntTransitThreshold >= 0.0);
    mWntTransitThreshold = wntTransitThreshold;
}

double SimpleWntCellCycleModel::GetWntLabelledThreshold() const
{
    return mWntLabelledThreshold;
}

void SimpleWntCellCycleModel::SetWntLabelledThreshold(double wntLabelledThreshold)
{
//    assert(wntLabelledThreshold <= 1.0);
//    assert(wntLabelledThreshold >= 0.0);
    mWntLabelledThreshold = wntLabelledThreshold;
}



void SimpleWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<UseCellProliferativeTypeDependentG1Duration>" << mUseCellProliferativeTypeDependentG1Duration << "</UseCellProliferativeTypeDependentG1Duration>\n";
    *rParamsFile << "\t\t\t<WntStemThreshold>" << mWntStemThreshold << "</WntStemThreshold>\n";
    *rParamsFile << "\t\t\t<WntTransitThreshold>" << mWntTransitThreshold << "</WntTransitThreshold>\n";
    *rParamsFile << "\t\t\t<WntLabelledThreshold>" << mWntLabelledThreshold << "</WntLabelledThreshold>\n";

    // Call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleWntCellCycleModel)
