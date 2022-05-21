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
#include "fold/bcl_fold_mutate_protein_model_sse_swap_multimer.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSESwapMultimer::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelSSESwapMultimer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from scheme
    //! @param BEND whether to bend the SSE after swapping to match the original phi/psi angles
    //! @param SCHEME scheme
    MutateProteinModelSSESwapMultimer::MutateProteinModelSSESwapMultimer
    (
      const bool BEND,
      const std::string &SCHEME
    ) :
      m_Bend( BEND),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSESwapMultimer
    MutateProteinModelSSESwapMultimer *MutateProteinModelSSESwapMultimer::Clone() const
    {
      return new MutateProteinModelSSESwapMultimer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSESwapMultimer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSESwapMultimer::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // cast a pointer to the multiplier data
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
      );

      // if the pointer is null
      if( !sp_multiplier.IsDefined())
      {
        // warn the user
        BCL_MessageStd
        (
          "Protein model multiplier not defined, unable to swap SSE with another subunit"
        );

        // return the normal score
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // generate the multimer
      const assemble::ProteinModel multimer_model( sp_multiplier->operator ()( PROTEIN_MODEL));

      // get the number of chains per subunit
      const size_t nr_chains( PROTEIN_MODEL.GetNumberOfChains());

      // get the number of subunits
      const size_t nr_subunits( multimer_model.GetNumberOfChains() / nr_chains);

      // get a random sse
      const util::SiPtrVector< const assemble::SSE> sses( PROTEIN_MODEL.GetSSEs());
      util::ShPtr< assemble::SSE> random_sse
      (
        ( *random::GetGlobalRandom().Iterator( sses.Begin(), sses.End(), sses.GetSize()))->Clone()
      );

      // get the chain id containing the sse to swap, random if it is before or after
      const util::ShPtr< assemble::Chain> swap_chain
      (
        random::GetGlobalRandom().Boolean() ?
          multimer_model.GetChain( char( random_sse->GetChainID() + nr_chains)) :
          multimer_model.GetChain( char( random_sse->GetChainID() + nr_chains * ( nr_subunits - 1)))
      );

      // get the sses in the chain
      util::SiPtrVector< const assemble::SSE> chain_sses( swap_chain->GetSSEs( random_sse->GetType()));

      // iterate through the sses
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( chain_sses.Begin()),
          sse_itr_end( chain_sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // if the SSE matches the random SSE
        if
        (
          ( *sse_itr)->GetFirstAA()->GetPdbID() == random_sse->GetFirstAA()->GetPdbID() &&
          ( *sse_itr)->GetLastAA()->GetPdbID() == random_sse->GetLastAA()->GetPdbID()
        )
        {
          // if the SSE should be bent
          if( m_Bend)
          {
            // fit the SSE
            random_sse->FitToSSE( **sse_itr);
          }
          else
          // don't bend
          {
            // move the random SSE to that position
            random_sse->Transform( math::Inverse( random_sse->GetOrientation()));
            random_sse->Transform( ( *sse_itr)->GetOrientation());
          }
          break;
        }
      }

      // make a copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // replace the random SSE
      new_model->Replace( random_sse);

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelSSESwapMultimer::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Bend, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelSSESwapMultimer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Bend, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace fold
  
} // namespace bcl
