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

#ifndef BCL_OPTI_EVOLUTION_OPERATION_H_
#define BCL_OPTI_EVOLUTION_OPERATION_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_evolution.h"
#include "bcl_opti_population.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EvolutionOperation
    //! @brief Interface for classes that determine when to terminate approximation
    //! @details ApproximatorModularBase derived classes will use EvolutionOperation derived classes to determine when
    //! to terminate their approximation. The criterion classes will be updated at each iteration and indicate when
    //! specific termination criteria are met
    //!
    //! @tparam t_MemberType is the type of the approximation argument
    //! @tparam t_FitnessType is the type of the approximation result
    //!
    //! @remarks example unnecessary
    //! @author geanesar
    //! @date Sep 16, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_MemberType, typename t_FitnessType>
    class EvolutionOperation :
      public util::FunctionInterfaceSerializable
      <
        util::SiPtrVector< const PopulationMember< t_MemberType, t_FitnessType> >,
        util::ShPtrVector< PopulationMember< t_MemberType, t_FitnessType> >
      >
    {

    public:

      //! @brief copy the class
      //! @return a pointer to a new instance of this class
      virtual EvolutionOperation *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief get the input size of the operation
      //! @return the number of members that must be in the SiPtrVector passed to operator()
      virtual size_t GetInputSize() const = 0;

    }; // template class EvolutionOperation< t_MemberType, t_FitnessType>
  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_EVOLUTION_OPERATION_H_
