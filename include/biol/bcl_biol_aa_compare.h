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

#ifndef BCL_BIOL_AA_COMPARE_H_
#define BCL_BIOL_AA_COMPARE_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_binary_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AACompareBySeqID
    //! @brief function class for comparison of a given amino acid by using sequence id
    //!
    //! @see @link example_biol_aa_compare.cpp @endlink
    //! @author karakam, alexanns
    //! @date Jan 22, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AACompareBySeqID
    {
    public:

    //////////
    // data //
    //////////

      const int m_SeqID; //!< sequence id

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a seq id
      //! @param SEQ_ID Sequence id
      AACompareBySeqID( const int SEQ_ID);

    ///////////////
    // operators //
    ///////////////

      //! @brief compares and returns whether the seq id of the given amino acid matched the stored one
      //! @param AMINO_ACID amino acoid to be compared
      //! @return returns whether the seq id of the given amino acid matched the stored one
      bool operator()( const AABase &AMINO_ACID) const;

      //! @brief compares and returns whether the seq id of the given amino acid matched the stored one
      //! @param SP_AMINO_ACID amino acoid to be compared
      //! @return returns whether the seq id of the given amino acid matched the stored one
      bool operator()( const util::PtrInterface< AABase> &SP_AMINO_ACID) const;

      //! @brief compares and returns whether the seq id of the given amino acid matched the stored one
      //! @param SP_AMINO_ACID amino acoid to be compared
      //! @return returns whether the seq id of the given amino acid matched the stored one
      bool operator()( const util::PtrInterface< const AABase> &SP_AMINO_ACID) const;

    }; // class AACompareBySeqID

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AACompareDataPtr
    //! @brief Compare two util::PtrInterface< const AABase> on whether the pointers are equal
    //!
    //! @see @link example_biol_aa_compare.cpp @endlink
    //! @author weinerbe, alexanns
    //! @date May 9, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AACompareDataPtr :
      public util::BinaryFunctionInterface< AABase, AABase, bool>
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new class_name
      AACompareDataPtr *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief () operator for returning whether two AABase are identical
      //! @param AA_BASE_LHS first AABase
      //! @param AA_BASE_RHS first AABase
      //! @return the whether two AABase are identical
      bool operator()( const AABase &AA_BASE_LHS, const AABase &AA_BASE_RHS) const;

      //! @brief () operator for returning whether two AABase are identical
      //! @param SP_AA_BASE_LHS first AABase
      //! @param SP_AA_BASE_RHS first AABase
      //! @return the whether two AABase are identical
      bool operator()
      (
        const util::PtrInterface< const AABase> &SP_AA_BASE_LHS,
        const util::PtrInterface< const AABase> &SP_AA_BASE_RHS
      ) const;

      //! @brief () operator for returning whether two AABase are identical
      //! @param SP_AA_BASE_LHS first AABase
      //! @param SP_AA_BASE_RHS first AABase
      //! @return the whether two AABase are identical
      bool operator()
      (
        const util::PtrInterface< AABase> &SP_AA_BASE_LHS,
        const util::PtrInterface< AABase> &SP_AA_BASE_RHS
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class AACompareDataPtr

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AACompareData
    //! @brief Compare two AABase on whether the stored data is identical
    //!
    //! @see @link example_biol_aa_compare.cpp @endlink
    //! @author weinerbe, alexanns
    //! @date May 9, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AACompareData :
      public util::BinaryFunctionInterface< AABase, AABase, bool>
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AACompareData()
      {
      }

      //! @brief Clone function
      //! @return pointer to new class_name
      AACompareData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief () operator for returning whether two AABase are identical
      //! @param AA_BASE_LHS first AABase
      //! @param AA_BASE_RHS first AABase
      //! @return the whether two AABase are identical
      bool operator()( const AABase &AA_BASE_LHS, const AABase &AA_BASE_RHS) const;

      //! @brief () operator for returning whether two AABase are identical
      //! @param SP_AA_BASE_LHS first AABase
      //! @param SP_AA_BASE_RHS first AABase
      //! @return the whether two AABase are identical
      bool operator()
      (
        const util::PtrInterface< const AABase> &SP_AA_BASE_LHS,
        const util::PtrInterface< const AABase> &SP_AA_BASE_RHS
      ) const;

      //! @brief () operator for returning whether two AABase are identical
      //! @param SP_AA_BASE_LHS first AABase
      //! @param SP_AA_BASE_RHS first AABase
      //! @return the whether two AABase are identical
      bool operator()
      (
        const util::PtrInterface< AABase> &SP_AA_BASE_LHS,
        const util::PtrInterface< AABase> &SP_AA_BASE_RHS
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class AACompareData

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AALessThanSeqID
    //! @brief Compare two AABase if one is less than the other based on chain id and seq id
    //!
    //! @see @link example_biol_aa_compare.cpp @endlink
    //! @author woetzen
    //! @date May 13, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AALessThanSeqID
    {
    public:

    ///////////////
    // operators //
    ///////////////

      //! @brief () operator for returning whether lhs amino acid is less than rhs amino acids by chain and seq id
      //! @param AA_BASE_LHS first AABase
      //! @param AA_BASE_RHS first AABase
      //! @return the whether one aa is less than the other aa
      bool operator()
      (
        const AABase &AA_BASE_LHS,
        const AABase &AA_BASE_RHS
      ) const;

      //! @brief () operator for returning whether lhs amino acid is less than rhs amino acids by chain and seq id
      //! @param SP_AA_BASE_LHS first AABase
      //! @param SP_AA_BASE_RHS first AABase
      //! @return the whether one aa is less than the other aa
      bool operator()
      (
        const util::PtrInterface< const AABase> &SP_AA_BASE_LHS,
        const util::PtrInterface< const AABase> &SP_AA_BASE_RHS
      ) const;

      //! @brief () operator for returning whether lhs amino acid is less than rhs amino acids by chain and seq id
      //! @param SP_AA_BASE_LHS first AABase
      //! @param SP_AA_BASE_RHS first AABase
      //! @return the whether one aa is less than the other aa
      bool operator()
      (
        const util::PtrInterface< AABase> &SP_AA_BASE_LHS,
        const util::PtrInterface< AABase> &SP_AA_BASE_RHS
      ) const;

    }; // class AALessThanSeqID

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_AA_COMPARE_H_ 
