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
#include "mc/bcl_mc_mutate_loop_add_resize.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_protein_geometry.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateLoopAddResize::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateLoopAddResize())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param LOOP_LIBRARY_FILENAME path of the loop template library
    //! @param MIN_SIZES minimum lengths of helices and strands after resize
    MutateLoopAddResize::MutateLoopAddResize
    (
      const std::string &LOOP_LIBRARY_FILENAME,
      const storage::VectorND< 2, size_t> MIN_SIZES
    ) :
      m_LoopLocator(),
      m_MinSizes( MIN_SIZES),
      m_LoopLibraryFilename( LOOP_LIBRARY_FILENAME)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief copy constructor
    //! @return pointer to a new MutateLoopAddResize
    MutateLoopAddResize *MutateLoopAddResize::Clone() const
    {
      return new MutateLoopAddResize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateLoopAddResize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default path of the loop template library
    //! @return the default path of the loop template library
    const std::string &MutateLoopAddResize::GetDefaultLibraryPath()
    {
      static std::string s_default_path( "histogram/loop_library");
      return s_default_path;
    }

    //! @brief returns the default minimum lengths of helices and strands after resize
    //! @return the default minimum lengths of helices and strands after resize
    storage::VectorND< 2, size_t> MutateLoopAddResize::GetDefaultMinSizes()
    {
      static storage::VectorND< 2, size_t> s_min_sizes( 5, 3);
      return s_min_sizes;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopAddResize::GetAlias() const
    {
      static const std::string s_name( "MutateLoopAddResize");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopAddResize::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Constructs a missing loop region in a protein model using a template library while randomly resizing the anchor SSEs."
      );
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
    bool MutateLoopAddResize::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_LoopLibrary = fold::LoopLibrary::CreateLoopLibrary( m_LoopLibraryFilename);
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief applies a mutation to the given protein model
    //! @detail this mutate randomly selects a missing loop in a protein model. the anchor SSEs of the selected loop are
    //! randomly resized and the resulting loop is constructed using a suitable conformation from a template library.
    //! @param MODEL protein model to which to apply the mutation
    //! @return result mutating the given protein model
    math::MutateResult< assemble::ProteinModel> MutateLoopAddResize::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // model resulting from applying the mutate
      util::ShPtr< assemble::ProteinModel> result_model( MODEL.Clone());

      // find missing loop regions in the protein model
      const util::ShPtrVector< fold::LoopParameters> loops( m_LoopLocator.Locate( MODEL));

      // if there are no missing loops in the protein model return the original model
      if( loops.IsEmpty())
      {
        BCL_MessageVrb( "No missing loops in the protein model. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // randomly select one of the missing loops for construction
      const size_t loop_index( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, loops.GetSize() - 1)));
      const fold::LoopParameters &loop( *loops( loop_index));

      // get the anchor SSEs of this loop
      const storage::VectorND< 2, util::ShPtr< assemble::SSE> > anchors( GetAnchors( MODEL, loop));

      // randomly resize the anchor SSEs
      const storage::VectorND< 2, util::ShPtr< assemble::SSE> > resized_sses( ResizeSSEs( anchors));
      const util::ShPtr< assemble::SSE> &sp_n_anchor( resized_sses.First());
      const util::ShPtr< assemble::SSE> &sp_c_anchor( resized_sses.Second());

      // compute the new loop for the resized anchors
      const util::ShPtr< fold::LoopParameters> sp_new_loop
      (
        fold::LoopParameters::Create( *sp_n_anchor->GetLastMember(), *sp_c_anchor->GetFirstMember())
      );

      // replace the anchor SSEs with the resized ones
      result_model->ReplaceResize( sp_n_anchor);
      result_model->ReplaceResize( sp_c_anchor);

      // find a suitable loop template in the library
      const util::ShPtrVector< fold::LoopParameters> templates( m_LoopLibrary->FindTemplates( loop));

      // loop region can't be constructed if there is no suitable template
      if( templates.IsEmpty())
      {
        // return start model if no suitable templates were found
        BCL_MessageVrb( "No suitable loop template found. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // select a random template from the found templates
      const util::ShPtrVector< fold::LoopParameters>::const_iterator selected_template_it
      (
        random::GetGlobalRandom().Iterator( templates.Begin(), templates.End(), templates.GetSize())
      );
      const fold::LoopParameters &selected_template( **selected_template_it);

      // fit the loop to the selected template and add it to the model
      result_model = fold::ProteinGeometry::FitToTemplate( MODEL, loop, selected_template);

      // return the result
      return math::MutateResult< assemble::ProteinModel>( result_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the anchor SSEs of a given loop
    //! @param MODEL protein model in which to find the anchor SSEs
    //! @param LOOP loop for which to return the anchor SSEs
    //! @return anchor SSEs of the given loop
    storage::VectorND< 2, util::ShPtr< assemble::SSE> > MutateLoopAddResize::GetAnchors
    (
      const assemble::ProteinModel &MODEL,
      const fold::LoopParameters &LOOP
    )
    {
      // get sequence and chain IDs of the anchor residues
      const storage::Vector< int> &anchor_ids( LOOP.GetAnchors());
      const char &chain_id( LOOP.GetChainID());

      // find the corresponding SSEs
      const assemble::Chain &chain( *MODEL.GetChain( chain_id));
      const util::SiPtrVector< const biol::AABase> residues( chain.GetSequence()->GetMembers());
      const biol::AABase &n_res( *residues( anchor_ids( 0) - 1));
      const biol::AABase &c_res( *residues( anchor_ids( 1) - 1));
      const util::ShPtr< assemble::SSE> sp_n_anchor( new assemble::SSE( *MODEL.GetSSE( n_res)));
      const util::ShPtr< assemble::SSE> sp_c_anchor( new assemble::SSE( *MODEL.GetSSE( c_res)));
      const storage::VectorND< 2, util::ShPtr< assemble::SSE> > anchors( sp_n_anchor, sp_c_anchor);

      return anchors;
    }

    //! @brief randomly resizes the given SSEs
    //! @param SSES the SSEs that shall be resized
    //! @return resized SSEs
    storage::VectorND< 2, util::ShPtr< assemble::SSE> > MutateLoopAddResize::ResizeSSEs
    (
      const storage::VectorND< 2, util::ShPtr< assemble::SSE> > &SSES
    )
    {
      // determine the sequence length of the loop and the length of the anchors
      const assemble::SSE &n_anchor( *SSES.First());
      const assemble::SSE &c_anchor( *SSES.Second());

      // determine which anchors to resize: 0 for both, 1 for n-terminal, 2 for c-terminal
      const size_t resize_targets( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, 2)));

      // determine the resize ranges
      const size_t resize_n
      (
        n_anchor.GetSize() > 6 ? random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, 2)) : 0
      );
      const size_t resize_c
      (
        c_anchor.GetSize() > 6 ? random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, 2)) : 0
      );

      // resize the SSEs
      storage::VectorND< 2, util::ShPtr< assemble::SSE> > resized_sses;
      // if( resize_targets == 0 || resize_targets == 2) // resize n-terminal SSE
      if( false)
      {
        resized_sses.First() = ShrinkSSE( n_anchor, resize_n, biol::AASequenceFlexibility::e_CTerminal);
      }
      else
      {
        resized_sses.First() = SSES.First();
      }
      if( resize_targets == 1 || resize_targets == 2) // resize c-terminal SSE
      {
        resized_sses.Second() = ShrinkSSE( c_anchor, resize_c, biol::AASequenceFlexibility::e_NTerminal);
      }
      else
      {
        resized_sses.Second() = SSES.Second();
      }

      return resized_sses;
    }

    //! @brief shrinks the given SSE by the given length
    //! @param SSE SSE to be shrunk
    //! @param LENGTH length by which the SSE will get shrunk
    //! @param DIRECTION terminus at which to shrink the SSE
    //! @return the shrunk SSE
    util::ShPtr< assemble::SSE> MutateLoopAddResize::ShrinkSSE
    (
      const assemble::SSE &SSE, size_t LENGTH, const biol::AASequenceFlexibility::SequenceDirection &DIRECTION
    )
    {
      // check if the SSE is longer than the size reduction
      const size_t new_length( SSE.GetSize() - LENGTH);
      BCL_Assert( new_length > 0, "Operation would result in SSE of negative length.");

      // create a new shrunk SSE
      util::ShPtr< assemble::SSE> sp_shrunk_sse;
      if( DIRECTION == biol::AASequenceFlexibility::e_NTerminal)
      {
        sp_shrunk_sse = util::ShPtr< assemble::SSE>
          (
            new assemble::SSE( SSE.SubSequence( LENGTH, SSE.GetLength()), SSE.GetType())
          );
      }
      else if( DIRECTION == biol::AASequenceFlexibility::e_CTerminal)
      {
        sp_shrunk_sse = util::ShPtr< assemble::SSE>
          (
            new assemble::SSE( SSE.SubSequence( 0, new_length), SSE.GetType())
          );
      }
      else
      {
        BCL_MessageCrt( "Direction " + util::Format()( DIRECTION) + " is unsupported. Returning original SSE.");
        sp_shrunk_sse = util::CloneToShPtr( SSE);
      }

      return sp_shrunk_sse;
    }

  } // namespace mc
} // namespace bcl
