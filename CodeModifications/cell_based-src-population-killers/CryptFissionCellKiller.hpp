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

#ifndef CRYPTFISSIONCELLKILLER_HPP_
#define CRYPTFISSIONCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "RandomNumberGenerator.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>



template<unsigned DIM>
class CryptFissionCellKiller : public AbstractCellKiller<DIM>
{
private:

    /**
     * Probability of crypt fission in an hour, elevated for mutant crypts.
      */
     double mCryptFissionRate;
	
	/** 
	 * Count number of crypts, increment after crypt fission.
	 */	
	 int crypt_counter;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);

        // Make sure the random number generator is also archived
        SerializableSingleton<RandomNumberGenerator>* p_rng_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_rng_wrapper;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param CryptFissionRate
     */
    CryptFissionCellKiller(AbstractCellPopulation<DIM>* pCellPopulation, double CryptFissionRate, int acrypt_counter);

    /**
     * @return mCryptFissionRate.
     */
    double GetCryptFissionRate() const;

	/** 
	 * @return crypt_counter.
	 */
	int GetCryptCounter() const;

    /**
     * Loop over cells and start apoptosis randomly, based on the user-set
     * probability.
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /**
     * Overridden OutputCellKillerParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CryptFissionCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptFissionCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CryptFissionCellKiller<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
    double prob = t->GetCryptFissionRate();
    ar << prob;
	int count = t-> GetCryptCounter();
	ar << count;
}

/**
 * De-serialize constructor parameters and initialise a CryptFissionCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CryptFissionCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
    double prob;
    ar >> prob;
	int count;
	ar >> count;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptFissionCellKiller<DIM>(p_cell_population, prob, count);
}
}
} // namespace ...

#endif /*CRYPTFISSIONCELLKILLER_HPP_*/
