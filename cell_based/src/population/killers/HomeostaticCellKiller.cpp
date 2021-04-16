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

#include "HomeostaticCellKiller.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "MMRTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
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


template<unsigned DIM>
HomeostaticCellKiller<DIM>::HomeostaticCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
        : AbstractCellKiller<DIM>(pCellPopulation)
{
}



template<unsigned DIM>
void HomeostaticCellKiller<DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
   
	// Start with homeostatic death not induced by cell size.

	double death_prob_this_timestep = 0.0;
	double ProbabilityOfDeathPerCC = 0.0005;
	double CellCycleTime = 24.0; // current TA cell cycle time, implemented in AbstractPhaseBasedCellCycleModel

	
	if (pCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>() || pCell->GetMutationState()->IsType<MMRApcTwoHitCellMutationState>() || pCell->GetMutationState()->IsType<BetaCateninTwoHitCellMutationState>() || pCell->GetMutationState()->IsType<MMRBCTwoHitCellMutationState>())
	{
		ProbabilityOfDeathPerCC = 0.0; // no homeostatic apoptosis for APC-- and beta-catenin++ mutated cells
	}	

	
	 /*
     * We have a constant time step dt, such that there are an integer number (n = 1/dt)
     * of time steps per hour.
	 * This means that there are m = CCT/dt time steps per cell cycle.
     *
     * Let q = death_prob_per_cell_cycle and p = death_prob_this_timestep
     *
     * Probability of not dying in one cell cycle:
     * (1-q) = (1-p)^m = (1-p)^(CCT/dt).
     *
     * Rearranging for p:
     * p = 1 - (1-q)^(dt/CCT).
     */

	death_prob_this_timestep = 1.0 - pow((1.0 - ProbabilityOfDeathPerCC), (SimulationTime::Instance()->GetTimeStep()/CellCycleTime));



    if (!pCell->HasApoptosisBegun() &&
        RandomNumberGenerator::Instance()->ranf() < death_prob_this_timestep)
    {
        pCell->StartApoptosis(); // cell dies
    }

	// Now implement the size constraint

	double cell_volume = pCell->GetCellData()->GetItem("volume"); // track cell volume
	pCell->RemoveCellProperty<CellLabel>();
	double size_threshold;

	if(pCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>() || pCell->GetMutationState()->IsType<MMRApcTwoHitCellMutationState>() || pCell->GetMutationState()->IsType<BetaCateninTwoHitCellMutationState>() || pCell->GetMutationState()->IsType<MMRBCTwoHitCellMutationState>())
	{
		size_threshold = 0.3; // APC-- and beta-catenin++ mutated cells respond less to apoptotic signals
	}
	else
	{
		size_threshold = 0.43;
	}


	if (!pCell->HasApoptosisBegun() && cell_volume < size_threshold)
	{	
		boost::shared_ptr<AbstractCellProperty> p_label = pCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CellLabel>();
		pCell->AddCellProperty(p_label);
		pCell->StartApoptosis();
	}
	

}

template<unsigned DIM>
void HomeostaticCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        CheckAndLabelSingleCellForApoptosis(*cell_iter);
    }
}


template<unsigned DIM>
void HomeostaticCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
   
    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class HomeostaticCellKiller<1>;
template class HomeostaticCellKiller<2>;
template class HomeostaticCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HomeostaticCellKiller)
