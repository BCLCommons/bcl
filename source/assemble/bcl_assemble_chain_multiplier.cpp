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
#include "assemble/bcl_assemble_chain_multiplier.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_chain.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ChainMultiplier::s_Instance
    (
      GetObjectInstances().AddInstance( new ChainMultiplier())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ChainMultiplier::ChainMultiplier() :
      m_Transformer(),
      m_InitialChainID( 'A'),
      m_Transformation(),
      m_Sequence()
    {
    }

    //! @brief construct from an SSE transformer, chain ID, and sequence
    //! @param SSE_TRANSFORMER transformer to apply to each SSE
    //! @param INITIAL_CHAIN_ID initial chain id
    //! @param TRANSFORMATION transformation that is applied to each SSE
    //! @param SEQUENCE sequence used to generate new chain
    ChainMultiplier::ChainMultiplier
    (
      const util::ShPtr< util::FunctionInterface< SSE, util::ShPtr< SSE> > > &SSE_TRANSFORMER,
      const char INITIAL_CHAIN_ID,
      const util::ShPtr< math::TransformationMatrix3D> &TRANSFORMATION,
      const util::ShPtr< biol::AASequence> &SEQUENCE
    ) :
      m_Transformer( SSE_TRANSFORMER),
      m_InitialChainID( INITIAL_CHAIN_ID),
      m_Transformation( TRANSFORMATION),
      m_Sequence( SEQUENCE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ChainMultiplier
    ChainMultiplier *ChainMultiplier::Clone() const
    {
      return new ChainMultiplier( *this);
    }

    //! @brief hardcopy this multiplier
    //! @return harcopied multiplier
    util::ShPtr< ChainMultiplier> ChainMultiplier::HardCopy() const
    {
      // set members
      util::ShPtr< ChainMultiplier> hard_copy( new ChainMultiplier());
      hard_copy->m_Transformer = m_Transformer.HardCopy();
      hard_copy->m_InitialChainID = m_InitialChainID;
      hard_copy->m_Transformation = m_Transformation.HardCopy();
      hard_copy->m_Sequence = util::ShPtr< biol::AASequence>( m_Sequence->HardCopy());

      // end
      return hard_copy;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ChainMultiplier::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief returns chain after transforming all SSEs
    //! @param CHAIN chain to be replicated
    //! @return chain after transforming all SSEs
    util::ShPtr< Chain> ChainMultiplier::operator ()( const Chain &CHAIN) const
    {
      // create new chain with the new sequence
      util::ShPtr< Chain> new_chain( new Chain( m_Sequence));

      // iterate through the SSEs on the chain
      util::SiPtrVector< const SSE> sses( CHAIN.GetSSEs());
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // insert the transformed SSE into the chain
        new_chain->Insert( m_Transformer->operator ()( **sse_itr));
      }

      // end
      return new_chain;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ChainMultiplier::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Transformer, ISTREAM);
      io::Serialize::Read( m_InitialChainID, ISTREAM);
      io::Serialize::Read( m_Transformation, ISTREAM);
      io::Serialize::Read( m_Sequence, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ChainMultiplier::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Transformer, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_InitialChainID, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Transformation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Sequence, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  /////////////////////////////
  // ChainMultiplierLessThan //
  /////////////////////////////

    //! @brief return whether one chain multiplier is less than another
    //! @param CHAIN_MULTIPLIER_LHS first chain multiplier
    //! @param CHAIN_MULTIPLIER_RHS second chain multiplier
    //! @return whether one chain multiplier is less than another
    bool ChainMultiplierLessThan::operator()
    (
      const ChainMultiplier &CHAIN_MULTIPLIER_LHS,
      const ChainMultiplier &CHAIN_MULTIPLIER_RHS
    ) const
    {
      // check if initial chain id of LHS is less than chain id of RHS
      if( CHAIN_MULTIPLIER_LHS.GetInitialChainID() < CHAIN_MULTIPLIER_RHS.GetInitialChainID())
      {
        return true;
      }

      // check if initial chain id of LHS is equal to chain id of RHS and new chain id of LHS is less than chain id
      // of RHS
      if
      (
        CHAIN_MULTIPLIER_LHS.GetInitialChainID() == CHAIN_MULTIPLIER_RHS.GetInitialChainID() &&
        CHAIN_MULTIPLIER_LHS.GetNewChainID() < CHAIN_MULTIPLIER_RHS.GetNewChainID()
      )
      {
        return true;
      }

      // else
      return false;
    }

    //! @brief return whether one chain multiplier is less than another
    //! @param PTR_CHAIN_MULTIPLIER_LHS first chain multiplier
    //! @param PTR_CHAIN_MULTIPLIER_RHS second chain multiplier
    //! @return whether one chain multiplier is less than another
    bool ChainMultiplierLessThan::operator()
    (
      const util::PtrInterface< ChainMultiplier> &PTR_CHAIN_MULTIPLIER_LHS,
      const util::PtrInterface< ChainMultiplier> &PTR_CHAIN_MULTIPLIER_RHS
    ) const
    {
      return operator ()( *PTR_CHAIN_MULTIPLIER_LHS, *PTR_CHAIN_MULTIPLIER_RHS);
    }

  } // namespace assemble
} // namespace bcl
