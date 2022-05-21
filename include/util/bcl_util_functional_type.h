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

#ifndef BCL_UTIL_FUNCTIONAL_TYPE_H_
#define BCL_UTIL_FUNCTIONAL_TYPE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FunctionalType
    //! @brief Logical type of a functor-style class
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Jan 28, 2015
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API FunctionalType
    {
    public:
      //! enum for what the functional type of an object is
      enum Type
      {
        e_StaticMethod,                 //!< No user-definable parameters (e.g. Atom_Identity)
        e_ParameterizedMethod,          //!< User-definable parameters (e.g. DihedralBins(bin size=30))
        // Wrapper (aka decorator design pattern) an object that derives from the same interface as a user-defined parameter of the object
        //   Typically the wrapper "decorates" or in some way modifies the output from the internal object.
        //   Wrappers are further classified by the number of decorators/objects of the same type that they use, and
        //   whether or not they have further parameters
        e_StaticWrapper,               //!< Takes an object of the same interface type, no other parameters (e.g. Sum)
        e_ParameterizedWrapper,        //!< Takes an object of the same interface type and other parameters (e.g. 3DA)
        s_NumberTypes
      };

      //! @brief get the type name
      //! @param TYPE the type name
      //! @return string describing the type
      static const std::string &GetTypeName( const Type &TYPE);

      //! @brief get the type from the number of wrappers and whether there are other parameters
      //! @param N_WRAPPERS the number of wrappers
      //! @param OTHER_PARAMS whether there were other parameters
      //! @return string describing the type
      static Type GetType( const size_t &N_WRAPPERS, const bool &OTHER_PARAMS);

      //! @brief detect whether a given type has non-decorator parameters
      //! @param TYPE the type
      //! @return bool - true if the type is parameterized
      static bool TestIsParameterized( const Type &TYPE);

      //! @brief detect the number of decorator parameters (3 means any number other than 0-2), e.g. 2 for Add, 1 for Sum
      //! @param TYPE the type
      //! @return number of user-selected classes with the same interface
      static size_t GetWrapperSize( const Type &TYPE);

    };
  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_FUNCTIONAL_TYPE_H_

