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
#include "mc/bcl_mc_mutate_loop_remove.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_complete.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> MutateLoopRemove::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateLoopRemove())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateLoopRemove::MutateLoopRemove() :
      m_LoopLocator( true, false) // locate defined loops and ignore terminal loops
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new MutateLoopRemove
    MutateLoopRemove *MutateLoopRemove::Clone() const
    {
      return new MutateLoopRemove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateLoopRemove::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopRemove::GetAlias() const
    {
      static const std::string s_name( "MutateLoopRemove");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopRemove::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Removes a loop from a protein model.");

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool MutateLoopRemove::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_LoopLocator = fold::LocatorLoop();
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief applies a mutation to the given protein model
    //! @param MODEL protein model to which to apply the mutation
    //! @return result mutating the given protein model
    math::MutateResult< assemble::ProteinModel> MutateLoopRemove::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // clone the original model
      util::ShPtr< assemble::ProteinModel> result_model( MODEL.Clone());

      // find defined loop regions in the model
      const util::ShPtrVector< fold::LoopParameters> loops( m_LoopLocator.Locate( MODEL));

      // if there are now defined loops in the protein model return the original model
      if( loops.IsEmpty())
      {
        BCL_MessageVrb( "No defined loops in the protein model. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // randomly select one of the loop regions
      const size_t loop_index( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, loops.GetSize() - 1)));
      const fold::LoopParameters &loop( *loops( loop_index));

      // find the loop in the protein model
      const util::ShPtr< assemble::Chain> &sp_chain( result_model->GetChain( loop.GetChainID()));
      const util::SiPtrVector< const biol::AABase> res_tmp( sp_chain->GetAminoAcids());
      const biol::AABase &first_loop_res( *res_tmp( loop.GetAnchors()( 0)));
      const assemble::SSE &selected_loop( *result_model->GetSSE( first_loop_res));
      const util::SiPtrVector< const biol::AABase> residues( selected_loop.GetMembers());

      // create a copy of the selected loop with undefined coordinates
      util::ShPtrVector< biol::AABase> new_seq;
      for( auto res_it( residues.Begin()), res_it_end( residues.End()); res_it != res_it_end; ++res_it)
      {
        const util::ShPtr< biol::AABase> sp_new_res( new biol::AAComplete( ( **res_it).GetData()));
        new_seq.PushBack( sp_new_res);
      }
      const biol::AASequence new_sequence( new_seq, selected_loop.GetFirstAA()->GetChainID());
      const util::ShPtr< assemble::SSE> sp_new_loop( new assemble::SSE( new_sequence, selected_loop.GetType()));

      // replace the selected loop in the protein model with the new undefined one
      result_model->Replace( sp_new_loop);
      util::ShPtr< assemble::SSE> sp_loop( new assemble::SSE( selected_loop, biol::GetSSTypes().COIL));

      // return the result
      const math::MutateResult< assemble::ProteinModel> mutate_result( result_model, *this);
      return mutate_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace mc
} // namespace bcl
