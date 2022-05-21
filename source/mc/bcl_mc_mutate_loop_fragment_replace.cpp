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
#include "mc/bcl_mc_mutate_loop_fragment_replace.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_protein_geometry.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateLoopFragmentReplace::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateLoopFragmentReplace)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param LOOP_LIBRARY_FILENAME path of the loop template library
    MutateLoopFragmentReplace::MutateLoopFragmentReplace( const std::string &LOOP_LIBRARY_FILENAME) :
      m_LoopLocator(),
      m_LoopLibraryFilename( LOOP_LIBRARY_FILENAME)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief copy constructor
    //! @return pointer to a new MutateLoopFragmentReplace
    MutateLoopFragmentReplace *MutateLoopFragmentReplace::Clone() const
    {
      return new MutateLoopFragmentReplace( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateLoopFragmentReplace::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default path of the loop template library
    //! @return the default path of the loop template library
    const std::string &MutateLoopFragmentReplace::GetDefaultLibraryPath()
    {
      static std::string s_default_path( "histogram/loop_library");
      return s_default_path;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopFragmentReplace::GetAlias() const
    {
      static const std::string s_name( "MutateLoopFragmentReplace");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopFragmentReplace::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Replaces a fragment of a loop in a protein model.");
      serializer.AddInitializer
      (
        "loop library",
        "path to the loop template library",
        io::Serialization::GetAgentInputFilename( &m_LoopLibraryFilename)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool MutateLoopFragmentReplace::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_LoopLibrary = fold::LoopLibrary::CreateLoopLibrary( m_LoopLibraryFilename);
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief applies a mutation to the given protein model
    //! @param MODEL protein model to which to apply the mutation
    //! @return result mutating the given protein model
    math::MutateResult< assemble::ProteinModel> MutateLoopFragmentReplace::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // model resulting from applying the mutate
      util::ShPtr< assemble::ProteinModel> result_model( MODEL.Clone());

      // find loop regions in the protein model
      const util::ShPtrVector< fold::LoopParameters> loops( m_LoopLocator.Locate( MODEL));
      util::ShPtrVector< fold::LoopParameters> loops_filtered;
      for( auto loop_it( loops.Begin()), loop_it_end( loops.End()); loop_it != loop_it_end; ++loop_it)
      {
        const fold::LoopParameters &loop_params( **loop_it);
        const char chain_id( loop_params.GetChainID());
        const storage::Pair< int, int> res_ids( loop_params.GetAnchors()( 0) + 1, loop_params.GetAnchors()( 1) - 1);
        const size_t length( res_ids.Second() - res_ids.First());
        const biol::AASequence loop( MODEL.GetChain( chain_id)->GetSequence()->SubSequence( res_ids.First() - 1, length));
        if( IsPartiallyDefined( loop))
        {
          loops_filtered.PushBack( *loop_it);
        }
      }

      // randomly select one fragment for reconstruction
      if( loops_filtered.IsEmpty())
      {
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }
      const math::Range< size_t> range( 0, loops_filtered.GetSize() - 1);
      const size_t loop_index( random::GetGlobalRandom().SizeT( range));
      const fold::LoopParameters &selected_loop( *loops_filtered( loop_index));
      const char chain_id( selected_loop.GetChainID());
      const storage::Pair< int, int> res_ids( selected_loop.GetAnchors()( 0), selected_loop.GetAnchors()( 1));
      const biol::AASequence loop_sequence
      (
        MODEL.GetChain( chain_id)->GetSequence()->SubSequence( res_ids.First(), res_ids.Second() - res_ids.First())
      );
      const util::ShPtr< fold::LoopParameters> selected_fragment( SelectFragment( loop_sequence));
      if( !selected_fragment.IsDefined())
      {
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }
      const size_t frag_length( selected_fragment->GetSequenceDistance() + 2);
      const math::Range< size_t> frag_range( 1, frag_length);
      const size_t selected_frag_length( random::GetGlobalRandom().SizeT( frag_range));
      storage::Vector< int> rebuild_ids;
      rebuild_ids.PushBack( selected_fragment->GetAnchors()( 1) - selected_frag_length);
      rebuild_ids.PushBack( selected_fragment->GetAnchors()( 1) + 1);
      const fold::LoopParameters fragment_rebuild( rebuild_ids, selected_fragment->GetChainID());

      // find a suitable template for the selected fragment in the library
      const util::ShPtrVector< fold::LoopParameters> templates( m_LoopLibrary->FindTemplates( selected_frag_length));
      if( templates.IsEmpty())
      {
        BCL_MessageTop( "No suitable loop template found. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }
      // select a random template from the found templates
      const util::ShPtrVector< fold::LoopParameters>::const_iterator selected_template_it
      (
        random::GetGlobalRandom().Iterator( templates.Begin(), templates.End(), templates.GetSize())
      );
      const fold::LoopParameters &selected_template( **selected_template_it);

      // fit the loop to the selected template and add it to the model
      result_model = fold::ProteinGeometry::FitToTemplate( MODEL, fragment_rebuild, selected_template);

      // return the result
      return math::MutateResult< assemble::ProteinModel>( result_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Returns whether the provided sequence is partially defined
    //! @param SEQUENCE sequence to be evaluated
    //! @return true, if the provided sequence is partially defined
    bool MutateLoopFragmentReplace::IsPartiallyDefined( const biol::AASequence &SEQUENCE) const
    {
      for( auto res_it( SEQUENCE.Begin()), res_it_end( SEQUENCE.End()); res_it != res_it_end; ++res_it)
      {
        if( ( **res_it).HasDefinedCoordinates())
        {
          for( auto res_int_it( res_it + 1); res_int_it != res_it_end; ++res_int_it)
          {
            if( ( **res_int_it).HasDefinedCoordinates())
            {
              return true;
            }
          }
        }
      }
      return false;
    }

    //! @brief Randomly select a fragment for reconstruction
    //! @param SEQUENCE sequence to select from
    //! @return randomly selected fragment
    util::ShPtr< fold::LoopParameters> MutateLoopFragmentReplace::SelectFragment( const biol::AASequence &SEQUENCE) const
    {
      const size_t n_id( ( **SEQUENCE.Begin()).GetSeqID());
      for( auto res_it( SEQUENCE.Begin() + 1), res_it_end( SEQUENCE.End()); res_it != res_it_end; ++res_it)
      {
        const biol::AABase &current_res( **res_it);
        if( !current_res.HasDefinedCoordinates())
        {
          storage::Vector< int> anchors;
          anchors.PushBack( n_id);
          anchors.PushBack( current_res.GetSeqID() - 1);
          const util::ShPtr< fold::LoopParameters> sp_frag( new fold::LoopParameters( anchors, current_res.GetChainID()));
          return sp_frag;
        }
      }
      return util::ShPtr< fold::LoopParameters>();
    }

  } // namespace mc
} // namespace bcl
