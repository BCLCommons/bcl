// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_OPTI_EVOLUTION_POPULATION_BUILDER_H_
#define BCL_OPTI_EVOLUTION_POPULATION_BUILDER_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_evolution.h"
#include "bcl_opti_evolution_population.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EvolutionPopulationBuilder
    //! @brief Interface for classes that determine when to terminate approximation
    //! @details ApproximatorModularBase derived classes will use EvolutionPopulationBuilder derived classes to determine when
    //! to terminate their approximation. The criterion classes will be updated at each iteration and indicate when
    //! specific termination criteria are met
    //!
    //! @tparam t_MemberType is the type of the approximation argument
    //! @tparam t_FitnessType is the type of the approximation result
    //!
    //! @remarks example unnecessary
    //! @author geanesar
    //! @date Jan 05, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionPopulationBuilder :
      public util::SerializableInterface
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief copy the class
      //! @return a pointer to a new instance of this class
      virtual EvolutionPopulationBuilder *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief build a new population from old populations
      //! @param POPULATION where to put the new members
      //! @param APPROXIMATOR the approximator that called this function
      virtual void Build
      (
        EvolutionPopulation< t_MemberType, t_FitnessType> &POPULATION,
        ApproximatorEvolution< t_MemberType, t_FitnessType> &APPROXIMATOR
      ) = 0;

    }; // template class EvolutionPopulationBuilder< t_MemberType, t_FitnessType>
  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_EVOLUTION_POPULATION_BUILDER_H_
