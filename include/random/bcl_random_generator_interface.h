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

#ifndef BCL_RANDOM_GENERATOR_INTERFACE_H_
#define BCL_RANDOM_GENERATOR_INTERFACE_H_

// include the namespace header
#include "bcl_random.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace random
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GeneratorInterface
    //! @brief is the interface for basic random number generators for uniform size_t distributions
    //! @details The SetSeed( SEED) and the SizeT() functions have to be overwritten in derived classes.
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Jul 9, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GeneratorInterface :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new GeneratorInterface
      virtual GeneratorInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief set the seed
      //! @param SEED to be used
      virtual void SetSeed( const uint64_t SEED) = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief generate random unsigned 64-bit integer
      //! @return integer uniformly distributed in range [0,2^64-1]
      virtual uint64_t Unsigned64BitInt() const = 0;

    }; // class GeneratorInterface

  } // namespace random
} // namespace bcl

#endif // BCL_RANDOM_GENERATOR_INTERFACE_H_ 
