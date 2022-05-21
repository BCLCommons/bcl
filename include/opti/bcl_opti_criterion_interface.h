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

#ifndef BCL_OPTI_CRITERION_INTERFACE_H_
#define BCL_OPTI_CRITERION_INTERFACE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_tracker.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CriterionInterface
    //! @brief Interface for classes that determine when to terminate approximation
    //! @details ApproximatorModularBase derived classes will use CriterionInterface derived classes to determine when
    //! to terminate their approximation. The criterion classes will be updated at each iteration and indicate when
    //! specific termination criteria are met
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @author fischea
    //! @date Dec 12, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionInterface :
      public util::SerializableInterface
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief Clone function
      //! @return pointer to a new CriterionInterface< t_ArgumentType, t_ResultType>
      virtual CriterionInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether any of the termination criteria are met
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, if any of the termination criteria are met
      virtual bool CriteriaMet( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const = 0;

    }; // template class CriterionInterface< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_CRITERION_INTERFACE_H_
