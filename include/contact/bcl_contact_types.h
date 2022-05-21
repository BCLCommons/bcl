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

#ifndef BCL_CONTACT_TYPES_H_
#define BCL_CONTACT_TYPES_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_contact_type_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Types
    //! @brief This is the enumerator class for representing contact types formed between residue pairs
    //! @details This enumerator class enumerates each contact type with a TypeData behind it. It also provides
    //! convenience functions to collect, reverse SSTypes and get valid distance ranges
    //!
    //! @see @link example_contact_types.cpp @endlink
    //! @author karakam
    //! @date 01.29.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Types :
      public util::Enumerate< TypeData, Types>
    {
      friend class util::Enumerate< TypeData, Types>;

    public:

    //////////
    // data //
    //////////

      // declare all secondary structure types
      const Type HELIX_HELIX;             //!< Helix-Helix SSType pairing
      const Type HELIX_SHEET;             //!< Helix-Sheet SSType pairing - when side chains of strand faces the helix
      const Type SHEET_HELIX;             //!< Sheet-Helix SSType pairing - when side chains of strand faces the helix
      const Type HELIX_STRAND;            //!< Helix-Strand SSType pairing - when hydrogens acceptors and donors of strand faces the helix
      const Type STRAND_HELIX;            //!< Strand-Helix SSType pairing - when hydrogens acceptors and donors of strand faces the helix
      const Type STRAND_STRAND;           //!< Strand-Strand SSType pairing - when hydrogen donor and acceptors of strands face each other
      const Type SHEET_SHEET;             //!< Sheet-Sheet SSType pairing   - when side chains of each strand face the other
      const Type UNDEFINED_HELIX_STRAND;  //!< undefined helix-strand pairing - twilight zone of side chain or hbond interactions
      const Type UNDEFINED_STRAND_HELIX;  //!< undefined strand-helix pairing - twilight zone of side chain or hbond interactions
      const Type UNDEFINED_STRAND_STRAND; //!< undefined strand-strand pairing - twilight zone of side chain or hbond interactions

      //! number of valid contact types
      enum { s_NumberValidTypes = 7};

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor that constructs all Types
      Types();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief cut off distance for unknown types;
      //! @return cut off distance for unknown types;
      static const double &GetUnknownResidueDistanceCutoff();

      //! @brief return distance range for unknown types
      //! @return distance range for unknown types
      static const math::Range< double> &GetUnknownDistanceRange();

      //! @brief return preferred distance range for unknown types
      //! @return preferred distance range for unknown types
      static const math::Range< double> &GetUnknownPreferredDistanceRange();

      //! @brief return tilt angle range for unknown types
      //! @return tilt angle range for unknown types
      static const math::Range< double> &GetUnknownTiltAngleRange();

      //! @brief return undefined pair length
      //! @return undefined pair length
      static const storage::Pair< size_t, size_t> &GetUndefinedLengthPair();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief returns the reverse pair for the contact type
      //! @param TYPE Type of interest
      //! @return the reverse pair for the contact type
      const Type &Reverse( const Type &TYPE) const;

      //! @brief Given two SSEGeometryInterfaces, merges their SSTypes to form a ContactType
      //! @param SSE_A first SSEGeometryInterface
      //! @param SSE_B second SSEGeometryInterface
      //! @return Contact type formed by SSE_A and SSE_B
      const Type &TypeFromSSTypes
      (
        const assemble::SSEGeometryInterface &SSE_A,
        const assemble::SSEGeometryInterface &SSE_B
      ) const;

      //! @brief returns a map with distance ranges for valid contact types
      //! @return a map with distance ranges for valid contact types
      const storage::Map< Type, math::Range< double> > &GetValidDistanceRanges() const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief generates a map with distance ranges for valid contact types
      //! @return a map with distance ranges for valid contact types
      storage::Map< Type, math::Range< double> > CollectValidDistanceRanges() const;

    }; // class Types

    //! @brief retrieves Types enumerator
    //! @return Types enumerator
    BCL_API const Types &GetTypes();

  } // namespace contact

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< contact::TypeData, contact::Types>;

  } // namespace util
} // namespace bcl

#endif //BCL_CONTACT_TYPES_H_
