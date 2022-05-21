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

#ifndef BCL_BIOL_AA_TYPES_H_
#define BCL_BIOL_AA_TYPES_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_type_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AATypes
    //! @brief This is the enumerator class amino acid types. Every enumerator has a distinct AATypeData behind it
    //! @details This enumerates allows representation of natural amino acid types as well as the unnatural amino acid
    //! types. Each amino acid type has a parent type which is itself for natural ones, and the corresponding natural
    //! type for the unnatural types.
    //!
    //! @see @link example_biol_aa_types.cpp @endlink
    //! @author woetzen, karakam
    //! @date 10.10.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AATypes :
      public util::Enumerate< AATypeData, AATypes>
    {
      friend class util::Enumerate< AATypeData, AATypes>;

    public:

    //////////
    // data //
    //////////

      // declare all amino acids
      AAType ALA; //!< Alanine
      AAType ARG; //!< Arginine
      AAType ASN; //!< Asparagine
      AAType ASP; //!< Aspartic acid
      AAType CYS; //!< Cysteine
      AAType GLN; //!< Glutamine
      AAType GLU; //!< Glutamic acid
      AAType GLY; //!< Glycine
      AAType HIS; //!< Histidine
      AAType ILE; //!< Isoleucine
      AAType LEU; //!< Leucine
      AAType LYS; //!< Lysine
      AAType MET; //!< Methionine
      AAType PHE; //!< Phenylalanine
      AAType PRO; //!< Proline
      AAType SER; //!< Serine
      AAType THR; //!< Threonine
      AAType TRP; //!< Tryptophan
      AAType TYR; //!< Tyrosine
      AAType VAL; //!< Valine

      AAType ASX; //!< Asparagine or Aspartic acid
      AAType GLX; //!< Glutamine or Glutamic acid
      AAType XXX; //!< arbitrary amino acid
      AAType UNK; //!< unknown amino acid as used in PDBs - do not mistake it for undefined AAType (UND)
      AAType GAP; //!< sequence gap

      AAType DAL; //!< D-ALANINE
      AAType ALS; //!< 2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
      AAType ACL; //!< DEOXY-CHLOROMETHYL-ARGININE
      AAType CME; //!< S,S-(2-HYDROXYETHYL)THIOCYSTEINE
      AAType CSE; //!< SELENOCYSTEINE
      AAType CSD; //!< 3-SULFINOALANINE
      AAType CSO; //!< S-HYDROXYCYSTEINE
      AAType CSW; //!< CYSTEINE-S-DIOXIDE
      AAType CYG; //!< 2-AMINO-4-(AMINO-3-OXO-PROPYLSULFANYLCARBONYL)-BUTYRIC ACID
      AAType OCS; //!< CYSTEINESULFONIC ACID
      AAType SC2; //!< N-ACETYL-L-CYSTEINE
      AAType CGU; //!< GAMMA-CARBOXY-GLUTAMIC ACID
      AAType PCA; //!< PYROGLUTAMIC ACID
      AAType KCX; //!< LYSINE NZ-CARBOXYLIC ACID
      AAType LLP; //!< 2-LYSINE(3-HYDROXY-2-METHYL-5-PHOSPHONOOXYMETHYL-PYRIDIN-4-YLMETHANE)
      AAType M3L; //!< N-TRIMETHYLLYSINE
      AAType MLY; //!< N-DIMETHYL-LYSINE
      AAType LYR; //!< N~6~-[(2Z,4E,6E,8E)-3,7-DIMETHYL-9-(2,6,6-TRIMETHYLCYCLOHEX-1-EN-1-YL)NONA-2,4,6,8-TETRAENYL]LYSINE
      AAType CXM; //!< N-CARBOXYMETHIONINE
      AAType FME; //!< N-FORMYLMETHIONINE
      AAType MSE; //!< Seleno Methionin
      AAType DPN; //!< D-PHENYLALANINE
      AAType S1H; //!< 1-HEXADECANOSULFONYL-O-L-SERINE
      AAType SAC; //!< N-ACETYL-SERINE
      AAType SEP; //!< PHOSPHOSERINE
      AAType TPO; //!< PHOSPHOTHREONINE
      AAType TYS; //!< O-SULFO-L-TYROSINE
      AAType R1A; //!< S-(2,2,5,5-tetramethyl-2,5-dihydro-1H-pyrrol-3-yl)methyl methanesulfonothioate
      AAType IAS; //!< BETA-L-ASPARTIC_ACID

      //! number of standard aminoacids
      enum { s_NumberStandardAATypes = 20};

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all AATypes
      AATypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief function to deduce AAType from three letter code of an amino acid
      //! @param THREE_LETTER_CODE three letter code descriptor for amino acid of interest
      //! @return AAType specified by given THREE_LETTER_CODE
      const AAType &AATypeFromThreeLetterCode( const std::string &THREE_LETTER_CODE) const;

      //! @brief function to deduce Parent AAType from three letter code of an amino acid
      //! this function is used to map unnatural amino acids back to their parents
      //! @param THREE_LETTER_CODE three letter code descriptor for amino acid of interest
      //! @return Parent for AAType specified by given THREE_LETTER_CODE
      const AAType &AATypeParentFromThreeLetterCode( const std::string &THREE_LETTER_CODE) const;

      //! @brief return if two amino acids have same parent
      //! this function returns true if two amino acids denoted by there three letter code have the same parent.
      //! MSE == MET would return true
      //! @param THREE_LETTER_CODE_LHS three letter code descriptor for amino acid of interest lhs
      //! @param THREE_LETTER_CODE_RHS three letter code descriptor for amino acid of interest rhs
      //! @return true if the two amino acids have same parent
      bool HaveSameParent( const std::string &THREE_LETTER_CODE_LHS, const std::string &THREE_LETTER_CODE_RHS) const;

      //! @brief function to deduce AAType from one letter code of an amino acid
      //! @param ONE_LETTER_CODE one letter code descriptor for amino acid of interest
      //! @return AAType specified by the given ONE_LETTER_CODE
      const AAType &AATypeFromOneLetterCode( const char ONE_LETTER_CODE) const;

      //! @brief gives the 20 natural amino acid types
      //! @return set of the 20 natural amino acid types
      const storage::Set< AAType> &GetNaturalAATypes() const;

    }; // class AATypes

    //! @brief construct on access function for all AATypes
    //! @return reference to only instance of AATypes enum
    BCL_API const AATypes &GetAATypes();

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< biol::AATypeData, biol::AATypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_BIOL_AA_TYPES_H_
