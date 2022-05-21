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

#ifndef BCL_BIOL_SS_TYPES_H_
#define BCL_BIOL_SS_TYPES_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_ss_type_data.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSTypes
    //! @brief This is the enumerator for different secondary structure types.
    //! @details Each instance of enum has specific SSTypeData behind it
    //!
    //! @see @link example_biol_ss_types.cpp @endlink
    //! @author karakam, woetzen
    //! @date 01.29.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSTypes :
      public util::Enumerate< SSTypeData, SSTypes>
    {
      friend class util::Enumerate< SSTypeData, SSTypes>;

    public:

    //////////
    // data //
    //////////

      // declare all secondary structure types
      const SSType HELIX;              //!< Helix right handed alpha; pdb helix class  1
      const SSType STRAND;             //!< Strand
      const SSType COIL;               //!< Coil
      const SSType e_HelixRightOmega;  //!< helix right handed omega; pdb helix class  2
      const SSType e_HelixRightPi;     //!< helix right handed pi   ; pdb helix class  3
      const SSType e_HelixRightGamma;  //!< helix right handed gamma; pdb helix class  4
      const SSType e_HelixRight310;    //!< helix right handed 3-10 ; pdb helix class  5
      const SSType e_HelixLeftAlpha;   //!< helix left handed alpha ; pdb helix class  6
      const SSType e_HelixLeftOmega;   //!< helix left handed omega ; pdb helix class  7
      const SSType e_HelixLeftGamma;   //!< helix left handed gamma ; pdb helix class  8
      const SSType e_Helix27Ribbon;    //!< helix 2 - 7 ribbon      ; pdb helix class  9
      const SSType e_HelixPolyProline; //!< helix poly proline      ; pdb helix class 10

    private:

      storage::Set< SSType> m_HelixTypes; //!< set of all helix types

    public:

      //! number of SSTypes between alpha helix and other helix types
      static const size_t s_AlphaOmegaHelixOffset = 2;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all SSTypes
      SSTypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the reduced (3-state) types
      const storage::Vector< SSType> &GetReducedTypes() const;

      //! @brief Get secondary structure type from provided ONE_LETTER_CODE
      //! @param ONE_LETTER_CODE one letter code for SSType
      //! @return secondary structure type from provided ONE_LETTER_CODE
      const SSType &SSTypeFromOneLetterCode( const char ONE_LETTER_CODE) const;

      //! @brief secondary structure from phi and psi angle
      //! @param PHI phi backbone angle [-pi,pi]
      //! @param PSI psi backbone angle [-pi,pi]
      //! @return most likely SSType - if non was found, COIL
      const SSType &SSTypeFromPhiPsi( const double PHI, const double PSI) const;

      //! @brief converts given pdb helix class into SSType
      //! @see http://www.wwpdb.org/documentation/format32/sect5.html#HELIX
      //! @param HELIX_CLASS as found in HelixClass entry in pdb HELIX line 1-10
      //! @return helix SSType for given class, e_Undefined if HELIX_CLASS is not in [1,10] range
      const SSType &SSTypeFromPDBHelixClass( const size_t HELIX_CLASS) const;

      //! @brief converts a SSType to a pdb helix class
      //! @param SSTYPE the sstype to convert to helix class
      //! @return pdb helix class [1-10]; if SSTYPE is not a valid helix class, undefined
      size_t PDBHelixClassFromSSType( const SSType &SSTYPE) const;

      //! @brief set of all helix types
      //! @return Set containing all helix types
      const storage::Set< SSType> &GetHelixTypes() const;

      //! @brief test if two sses are of similar type (alpha helix to 3-10 helix are considered similar, but strand would not be)
      //! @param SS_TYPE_LHS left hand side ss type
      //! @param SS_TYPE_RHS right hand side ss type
      //! @return true if the two types are considered similar
      bool AreSimilar( const SSType &SS_TYPE_LHS, const SSType &SS_TYPE_RHS) const;

    }; // class SSTypes

    //! @brief on access function for all SSTypes
    //! @return SSTypes enums
    BCL_API const SSTypes &GetSSTypes();

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< biol::SSTypeData, biol::SSTypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_BIOL_SS_TYPES_H_
