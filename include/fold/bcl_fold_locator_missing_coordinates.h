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

#ifndef BCL_FOLD_LOCATOR_MISSING_COORDINATES_H_
#define BCL_FOLD_LOCATOR_MISSING_COORDINATES_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorMissingCoordinates
    //!
    //! @see @link example_fold_locator_missing_coordinates.cpp @endlink
    //! @author fischea
    //! @date Mar 3, 2017
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API LocatorMissingCoordinates :
      public util::SerializableInterface
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! sequence span having missing backbone coordinates
      typedef storage::Triplet< int, int, char> Span;

    private:

      //! ignore terminal loops
      bool m_IgnoreTerms;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief construct from members
      //! @param IGNORE_TERMS ignore terminal loops
      LocatorMissingCoordinates( bool IGNORE_TERMS = true);

      //! @brief clone function
      //! @return pointer to a new LocatorMissingCoordinates
      LocatorMissingCoordinates *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief return regions with undefined backbone coordinates in the given protein model
      //! @param MODEL protein model for which to return regions with undefined backbone coordinates
      //! @return regions with undefined backbone coordinates
      storage::Vector< Span> Locate( const assemble::ProteinModel &MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief return regions with undefined backbone coordinates in the given sequence
      //! @param CHAIN sequence for which to return regions with undefined backbone coordinates
      //! @return regions with undefined backbone coordinates
      storage::Vector< LocatorMissingCoordinates::Span> Locate( const biol::AASequence &CHAIN) const;

    }; // class LocatorMissingCoordinates

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_LOCATOR_MISSING_COORDINATES_H_
