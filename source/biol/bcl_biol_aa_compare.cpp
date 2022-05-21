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
#include "biol/bcl_biol_aa_compare.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////
  // AACompareBySeqID //
  //////////////////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a seq id
    //! @param SEQ_ID Sequence id
    AACompareBySeqID::AACompareBySeqID( const int SEQ_ID) :
      m_SeqID( SEQ_ID)
    {
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief compares and returns whether the seq id of the given amino acid matched the stored one
    //! @param AMINO_ACID amino acoid to be compared
    //! @return returns whether the seq id of the given amino acid matched the stored one
    bool AACompareBySeqID::operator()( const AABase &AMINO_ACID) const
    {
      return AMINO_ACID.GetSeqID() == m_SeqID;
    }

    //! @brief compares and returns whether the seq id of the given amino acid matched the stored one
    //! @param SP_AMINO_ACID amino acoid to be compared
    //! @return returns whether the seq id of the given amino acid matched the stored one
    bool AACompareBySeqID::operator()( const util::PtrInterface< AABase> &SP_AMINO_ACID) const
    {
      return operator()( *SP_AMINO_ACID);
    }

    //! @brief compares and returns whether the seq id of the given amino acid matched the stored one
    //! @param SP_AMINO_ACID amino acoid to be compared
    //! @return returns whether the seq id of the given amino acid matched the stored one
    bool AACompareBySeqID::operator()( const util::PtrInterface< const AABase> &SP_AMINO_ACID) const
    {
      return operator()( *SP_AMINO_ACID);
    }

  //////////////////////
  // AACompareDataPtr //
  //////////////////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new class_name
    AACompareDataPtr *AACompareDataPtr::Clone() const
    {
      return new AACompareDataPtr( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AACompareDataPtr::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator for returning whether two AABase are identical
    //! @param AA_BASE_LHS first AABase
    //! @param AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareDataPtr::operator()( const AABase &AA_BASE_LHS, const AABase &AA_BASE_RHS) const
    {
      return AA_BASE_LHS.GetData() == AA_BASE_RHS.GetData();
    }

    //! @brief () operator for returning whether two AABase are identical
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareDataPtr::operator()
    (
      const util::PtrInterface< const AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< const AABase> &SP_AA_BASE_RHS
    ) const
    {
      return operator ()( *SP_AA_BASE_LHS, *SP_AA_BASE_RHS);
    }

    //! @brief () operator for returning whether two AABase are identical
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareDataPtr::operator()
    (
      const util::PtrInterface< AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< AABase> &SP_AA_BASE_RHS
    ) const
    {
      return operator ()( *SP_AA_BASE_LHS, *SP_AA_BASE_RHS);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AACompareDataPtr::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AACompareDataPtr::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  ///////////////////
  // AACompareData //
  ///////////////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new class_name
    AACompareData *AACompareData::Clone() const
    {
      return new AACompareData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AACompareData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator for returning whether two AABase are identical
    //! @param AA_BASE_LHS first AABase
    //! @param AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareData::operator()( const AABase &AA_BASE_LHS, const AABase &AA_BASE_RHS) const
    {
      return AA_BASE_LHS.GetChainID() == AA_BASE_RHS.GetChainID() &&
        AA_BASE_LHS.GetSeqID() == AA_BASE_RHS.GetSeqID() &&
        AA_BASE_LHS.GetType() == AA_BASE_RHS.GetType();
    }

    //! @brief () operator for returning whether two AABase are identical
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareData::operator()
    (
      const util::PtrInterface< const AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< const AABase> &SP_AA_BASE_RHS
    ) const
    {
      return operator ()( *SP_AA_BASE_LHS, *SP_AA_BASE_RHS);
    }

    //! @brief () operator for returning whether two AABase are identical
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether two AABase are identical
    bool AACompareData::operator()
    (
      const util::PtrInterface< AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< AABase> &SP_AA_BASE_RHS
    ) const
    {
      return operator ()( *SP_AA_BASE_LHS, *SP_AA_BASE_RHS);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AACompareData::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AACompareData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  /////////////////////
  // AALessThanSeqID //
  /////////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator for returning whether lhs amino acid is less than rhs amino acids by chain and seq id
    //! @param AA_BASE_LHS first AABase
    //! @param AA_BASE_RHS first AABase
    //! @return the whether one aa is less than the other aa
    bool AALessThanSeqID::operator()
    (
      const AABase &AA_BASE_LHS,
      const AABase &AA_BASE_RHS
    ) const
    {
      // compare chain id
      if( AA_BASE_LHS.GetChainID() == AA_BASE_RHS.GetChainID())
      {
        return AA_BASE_LHS.GetSeqID() < AA_BASE_RHS.GetSeqID();
      }

      return AA_BASE_LHS.GetChainID() < AA_BASE_RHS.GetChainID();
    }

    //! @brief () operator for returning whether lhs amino acid is less than rhs amino acids by chain and seq id
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether one aa is less than the other aa
    bool AALessThanSeqID::operator()
    (
      const util::PtrInterface< const AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< const AABase> &SP_AA_BASE_RHS
    ) const
    {
      // compare chain id
      if( SP_AA_BASE_LHS->GetChainID() == SP_AA_BASE_RHS->GetChainID())
      {
        return SP_AA_BASE_LHS->GetSeqID() < SP_AA_BASE_RHS->GetSeqID();
      }

      return SP_AA_BASE_LHS->GetChainID() < SP_AA_BASE_RHS->GetChainID();
    }

    //! @brief () operator for returning whether lhs amino acid is less than rhs amino acids by chain and seq id
    //! @param SP_AA_BASE_LHS first AABase
    //! @param SP_AA_BASE_RHS first AABase
    //! @return the whether one aa is less than the other aa
    bool AALessThanSeqID::operator()
    (
      const util::PtrInterface< AABase> &SP_AA_BASE_LHS,
      const util::PtrInterface< AABase> &SP_AA_BASE_RHS
    ) const
    {
      // compare chain id
      if( SP_AA_BASE_LHS->GetChainID() == SP_AA_BASE_RHS->GetChainID())
      {
        return SP_AA_BASE_LHS->GetSeqID() < SP_AA_BASE_RHS->GetSeqID();
      }

      return SP_AA_BASE_LHS->GetChainID() < SP_AA_BASE_RHS->GetChainID();
    }

  } // namespace biol
} // namespace bcl
