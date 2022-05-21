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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "math/bcl_math.h"
#include "opti/bcl_opti_improvement_type.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    //! @brief returns ImprovementType as string
    //! @param TYPE the improvement type
    //! @return the string for the improvement type
    const std::string &GetImprovementTypeName( const ImprovementType &TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "SmallerIsBetter",         //!< smaller values indicate improvement
        "LargerIsBetter",          //!< larger values indicate improvement
        "SmallerAbsIsBetter",      //!< smaller absolute values indicate improvement (typically error)
        "SmallerEqualIsBetter",    //!< smaller or equal values indicate improvement
        "LargerEqualIsBetter",     //!< larger or equal values indicate improvement
        GetStaticClassName< ImprovementType>()
      };

      return s_descriptors[ TYPE];
    }

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    template< typename t_DataType>
    bool TestDoesImprove( const t_DataType &TEST, const t_DataType &PRIOR, const ImprovementType &TYPE)
    {
      bool improved( false);
      switch( TYPE)
      {
        case e_SmallerIsBetter:      improved = TEST <  PRIOR; break;
        case e_LargerIsBetter:       improved = TEST >  PRIOR; break;
        case e_SmallerEqualIsBetter: improved = TEST <= PRIOR; break;
        case e_LargerEqualIsBetter:  improved = TEST >= PRIOR; break;
        case e_SmallerAbsIsBetter:   improved = math::Absolute( TEST) < math::Absolute( PRIOR); break;
        default: break;
      }
      return improved;
    }

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    bool DoesImprove( const double &TEST, const double &PRIOR, const ImprovementType &TYPE)
    {
      return TestDoesImprove( TEST, PRIOR, TYPE);
    }

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    bool DoesImprove( const float &TEST, const float &PRIOR, const ImprovementType &TYPE)
    {
      return TestDoesImprove( TEST, PRIOR, TYPE);
    }

    //! @brief test whether a particular value should be considered improved, given the improvement type
    //! @param TEST the value to test for improvement
    //! @param PRIOR the prior best value
    //! @param TYPE the improvement type
    //! @return true if TEST improves PRIOR based on TYPE
    bool DoesImprove( const int &TEST, const int &PRIOR, const ImprovementType &TYPE)
    {
      return TestDoesImprove( TEST, PRIOR, TYPE);
    }

  } // namespace opti
} // namespace bcl
