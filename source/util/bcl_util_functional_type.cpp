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
#include "util/bcl_util_functional_type.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the type name
    //! @param TYPE the type name
    //! @return string describing the type
    const std::string &FunctionalType::GetTypeName( const FunctionalType::Type &TYPE)
    {
      static const std::string s_names[ s_NumberTypes + 1] =
      {
        "Basic Implementations",        // No user-definable parameters (e.g. Atom_Identity)
        "Customizable Implementations", // User-definable parameters (e.g. DihedralBins(bin size=30))
        "Operations",                   // Takes an object of the same interface type, no other parameters (e.g. Sum)
        "Customizable Operations",      // Takes an object of the same interface type and other parameters (e.g. 3DA)
        GetStaticClassName< FunctionalType::Type>()
      };
      return s_names[ TYPE];
    }

    //! @brief get the type from the number of wrappers and whether there are other parameters
    //! @param N_WRAPPERS the number of wrappers
    //! @param OTHER_PARAMS whether there were other parameters
    //! @return string describing the type
    FunctionalType::Type FunctionalType::GetType( const size_t &N_WRAPPERS, const bool &OTHER_PARAMS)
    {
      return Type( ( N_WRAPPERS ? 2 : 0) + ( OTHER_PARAMS ? 1 : 0));
    }

    //! @brief detect whether a given type has non-decorator parameters
    //! @param TYPE the type
    //! @return bool - true if the type is parameterized
    bool FunctionalType::TestIsParameterized( const Type &TYPE)
    {
      return size_t( TYPE) & 1;
    }

    //! @brief detect the number of decorator parameters (3 means any number other than 0-2), e.g. 2 for Add, 1 for Sum
    //! @param TYPE the type
    //! @return number of user-selected classes with the same interface
    size_t FunctionalType::GetWrapperSize( const Type &TYPE)
    {
      return size_t( TYPE) / 2;
    }

  } // namespace util
} // namespace bcl
