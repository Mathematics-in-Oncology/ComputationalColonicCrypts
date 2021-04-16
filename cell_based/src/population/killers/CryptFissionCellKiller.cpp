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

#include "CryptFissionCellKiller.hpp"
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
#include "BCLOHApcOneHitCellMutationState.hpp"
#include "BCLOHMMRTwoHitCellMutationState.hpp"
#include "BCLOHMMR_ApcOneHitCellMutationState.hpp"
#include "BCLOHApcLOHCellMutationState.hpp"
#include "BCLOHMMR_ApcLOHCellMutationState.hpp"


template<unsigned DIM>
CryptFissionCellKiller<DIM>::CryptFissionCellKiller(AbstractCellPopulation<DIM>* pCellPopulation, double CryptFissionRate, int acrypt_counter)
        : AbstractCellKiller<DIM>(pCellPopulation),
          mCryptFissionRate(CryptFissionRate),
		  crypt_counter(acrypt_counter)
{
    if ((mCryptFissionRate<0) || (mCryptFissionRate>1))
    {
        EXCEPTION("Probability of crypt fission must be between zero and one");
    }
}

template<unsigned DIM>
double CryptFissionCellKiller<DIM>::GetCryptFissionRate() const
{
    return mCryptFissionRate;
}

template<unsigned DIM>
int CryptFissionCellKiller<DIM>::GetCryptCounter() const
{
    return crypt_counter;
}


template<unsigned DIM>
void CryptFissionCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{	

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	// initialize MPI rank
	
    /*
     * We assume a constant time step and that there are an integer number (n = 1/dt)
     * of time steps per hour. We also assume that this method is called every time step
     * and that the probabilities of fission at different times are independent.
     *
     * Let q=mCryptFissionRate and p="probability of fission in a given time step".
     *
     * Probability of no fission in an hour:
     * (1-q) = (1-p)^n = (1-p)^(1/dt).
     *
     * Rearranging for p:
     * p = 1 - (1-q)^dt.
     */

	double fission_prob_this_timestep = 0.0;
	double mmr_percentage = (double)this->mpCellPopulation->GetNumMMRCells()/(double)this->mpCellPopulation->GetNumRealCells();

	// We declare a crypt as MMR-deficient if over 80 percent of all cells show a double MMR-hit.
	// Crypt fission is set to be more likely for such crypts.

	if (mmr_percentage > 0.8) 
	{
		//Increase fission probability of MMR-mutant crypts.
    	fission_prob_this_timestep = 1.0 - pow((1.0 - 1.22*mCryptFissionRate), SimulationTime::Instance()->GetTimeStep());
	}
	else 
	{
		fission_prob_this_timestep = 1.0 - pow((1.0 - mCryptFissionRate), SimulationTime::Instance()->GetTimeStep());
	}

    if (RandomNumberGenerator::Instance()->ranf() < fission_prob_this_timestep) // crypt fission occurs
    { 	
		crypt_counter++;  // add (virtual) crypt
		mCryptFissionRate = (double) crypt_counter * mCryptFissionRate/((double)crypt_counter -1.0);  // new fission rate is no. of crypts * inital fission rate
		if (mmr_percentage > 0.8)
		{			
			cout << "Crypt fission of MMR-deficient crypt in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl;   // output time of mutated crypt fission event   
		}
		else
		{
			cout << "Crypt fission of normal crypt in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() <<  endl;   // output time of normal crypt fission event
		} 
    }

	/*
	// mice simulation

	fission_prob_this_timestep = 1.0 - pow((1.0 - 1.22*mCryptFissionRate), SimulationTime::Instance()->GetTimeStep());
	if (RandomNumberGenerator::Instance()->ranf() < fission_prob_this_timestep) // crypt fission occurs
    { 	
		crypt_counter++;  // add (virtual) crypt
		mCryptFissionRate = (double) crypt_counter * mCryptFissionRate/((double)crypt_counter -1.0);  // new fission rate is no. of crypts * inital fission rate		
		cout << "Crypt fission of MMR-deficient crypt in cluster " << rank << " at time " << SimulationTime::Instance()->GetTime() << endl;   // output time of mutated crypt fission event
	}*/
}

template<unsigned DIM>
void CryptFissionCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CryptFissionRate>" << mCryptFissionRate << "</CryptFissionRate>\n";

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class CryptFissionCellKiller<1>;
template class CryptFissionCellKiller<2>;
template class CryptFissionCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CryptFissionCellKiller)
