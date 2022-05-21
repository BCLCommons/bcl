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

#ifndef BCL_ASSEMBLE_AA_EXPOSURES_H_
#define BCL_ASSEMBLE_AA_EXPOSURES_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_aa_exposure_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAExposures
    //! @brief enumerates aa exposure interfaces
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_assemble_aa_exposures.cpp @endlink
    //! @author woetzen
    //! @date Feb 13, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAExposures :
      public util::Enumerate< util::ShPtr< AAExposureInterface>, AAExposures>
    {
      friend class util::Enumerate< util::ShPtr< AAExposureInterface>, AAExposures>;
    public:

    //////////
    // data //
    //////////

      AAExposure e_NeighborCount;
      AAExposure e_NeighborVector;
      AAExposure e_AASasaOls;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      AAExposures();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    }; // class AAExposures

    //! @brief construct on access function for all AAExposures
    //! @return reference to only instances of AAExposures
    BCL_API
    AAExposures &GetAAExposures();

  } // namespace assemble

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< assemble::AAExposureInterface>, assemble::AAExposures>;

  } // namespace util
} // namespace bcl

#endif // BCL_ASSEMBLE_AA_EXPOSURES_H_ 
