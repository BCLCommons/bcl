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

#ifndef BCL_OPTI_IMPROVEMENT_TYPE_H_
#define BCL_OPTI_IMPROVEMENT_TYPE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

// includes from bcl - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @file bcl_opti_improvement_type.h
    //! @brief enum wrapper indicating how to test whether an optimization criterion improved
    //!
    //! @see @link example_opti_improvement_type.cpp @endlink
    //! @author mendenjl
    //! @date Aug 22, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! enumerator for change types that indicate an improvement
    enum ImprovementType
    {
      e_SmallerIsBetter,         //!< smaller values indicate improvement
      e_LargerIsBetter,          //!< larger values indicate improvement
      e_SmallerAbsIsBetter,      //!< smaller absolute values indicate improvement (typically error)
      e_SmallerEqualIsBetter,    //!< smaller or equal values indicate improvement
      e_LargerEqualIsBetter,     //!< larger or equal values indicate improvement
      s_NumberImprovementTypes
    };

    //! @brief returns ImprovementType as string
    //! @param TYPE the improvement type
    //! @return the string for the improvement type
    const std::string &GetImprovementTypeName( const ImprovementType &TYPE);

    //! SSEInfoTypeEnum simplifies the usage of the SSEInfoType enum of this class
    typedef util::WrapperEnum< ImprovementType, &GetImprovementTypeName, s_NumberImprovementTypes> ImprovementTypeEnum;

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    bool DoesImprove( const double &TEST, const double &PRIOR, const ImprovementType &TYPE);

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    bool DoesImprove( const float &TEST, const float &PRIOR, const ImprovementType &TYPE);

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    bool DoesImprove( const int &TEST, const int &PRIOR, const ImprovementType &TYPE);

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_IMPROVEMENT_TYPE_H_
