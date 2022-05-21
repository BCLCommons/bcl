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
#include "fold/bcl_fold_mutate_protein_model_sse_resize.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSEResize::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModelSSEResize())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSEResize::MutateProteinModelSSEResize() :
      m_SSELocator(),
      m_ExtendProbability( 0.5),
      m_LengthChangeRange(),
      m_Side(),
      m_RecenterAfterResize(),
      m_MinSSESizes(),
      m_Scheme( GetStaticClassName< MutateProteinModelSSEResize>())
    {
    }

    //! @brief constructor from a locator, an extend/shrink probability, length increment and a boolean flag
    //! @param SSE_LOCATOR locator that decides which sse in the proteinmodel to mutate
    //! @param EXTEND_PROBABILITY probability for extending (1.0 - EXTEND_PROBABILITY for shrinking)
    //! @param LENGTH_CHANGE_RANGE range of number of residues that are to be added or removed in one mutate to one end
    //! @param SEQUENCE_SIDE side of sequence to modify
    //! @param RECENTER_AFTER_RESIZE boolean to whether recenter the sse after resize to its original center
    //! @param MIN_SSE_SIZES map of minimum SSE sizes to be allowed when shrinking
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEResize::MutateProteinModelSSEResize
    (
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &SSE_LOCATOR,
      const double &EXTEND_PROBABILITY,
      const math::Range< size_t> &LENGTH_CHANGE_RANGE,
      const biol::AASequenceFlexibility::SequenceDirection &SEQUENCE_SIDE,
      const bool RECENTER_AFTER_RESIZE,
      const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES,
      const std::string &SCHEME
    ) :
      m_SSELocator( *SSE_LOCATOR),
      m_ExtendProbability( EXTEND_PROBABILITY),
      m_LengthChangeRange( LENGTH_CHANGE_RANGE),
      m_Side( SEQUENCE_SIDE),
      m_RecenterAfterResize( RECENTER_AFTER_RESIZE),
      m_MinSSESizes( MIN_SSE_SIZES),
      m_Scheme( SCHEME)
    {
      static const math::Range< double> s_default_range( 0.0, 1.0);
      // check that probability is within range
      BCL_Assert
      (
        s_default_range.IsWithin( m_ExtendProbability),
        "The given probability should be within range: " + util::Format()( s_default_range.GetString()) + " " +
        util::Format()( EXTEND_PROBABILITY)
      );
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSEResize
    MutateProteinModelSSEResize *MutateProteinModelSSEResize::Clone() const
    {
      return new MutateProteinModelSSEResize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelSSEResize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return extend probability
    //! @return extend probability
    const double &MutateProteinModelSSEResize::GetExtendProbability() const
    {
      return m_ExtendProbability;
    }

    //! @brief return length change range
    //! @return length change range
    const math::Range< size_t> &MutateProteinModelSSEResize::GetLengthChangeRange() const
    {
      return m_LengthChangeRange;
    }

    //! @brief return which ends are changed
    //! @return sequence direction
    biol::AASequenceFlexibility::SequenceDirection MutateProteinModelSSEResize::GetSide() const
    {
      return m_Side;
    }

    //! @brief returns min sse sizes
    //! @return min sse sizes
    const storage::Map< biol::SSType, size_t> &MutateProteinModelSSEResize::GetMinSSESizes() const
    {
      return m_MinSSESizes;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSEResize::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model
      const math::MutateResult< assemble::ProteinModel> empty_result( util::ShPtr< assemble::ProteinModel>(), *this);

      // locate an sse and make copy
      const util::SiPtr< const assemble::SSE> located_sse( m_SSELocator->Locate( PROTEIN_MODEL));

      // if sse cannot be located or the located sse is not helix or strand
      if( !located_sse.IsDefined() || located_sse->GetType() > biol::GetSSTypes().STRAND)
      {
        // return empty result
        return empty_result;
      }

      // report selected sse to be extended/shrunk
      BCL_MessageVrb( "selected sse to be extended/shrunk: " + located_sse->GetIdentification());

      // determine whether to extend or shrink
      const bool extend( random::GetGlobalRandom().Double() <= m_ExtendProbability);

      // determine step size
      const size_t step_size( random::GetGlobalRandom().Random< size_t>( m_LengthChangeRange));

      // make a ShPtr to an SSE
      util::ShPtr< assemble::SSE> new_sse;

      // if shrinking is picked
      if( !extend)
      {
        // find the min sse size
        const storage::Map< biol::SSType, size_t>::const_iterator itr( m_MinSSESizes.Find( located_sse->GetType()));
        const size_t min_sse_size( itr == m_MinSSESizes.End() ? 0 : itr->second);

        // is sse long enough
        if( m_Side == biol::AASequenceFlexibility::e_Bidirectional ? 2 * step_size > located_sse->GetSize() : step_size > located_sse->GetSize())
        {
          return empty_result;
        }

        // calculate the size of SSE after the shrinking
        const size_t length_after_shrink
        (
          m_Side == biol::AASequenceFlexibility::e_Bidirectional ? located_sse->GetSize() - 2 * step_size : located_sse->GetSize() - step_size
        );

        // remaining sse long enough?
        if( length_after_shrink < min_sse_size)
        {
          return empty_result;
        }

        // where to shrink
        switch( m_Side)
        {
          // if both ends are shrunk
          case biol::AASequenceFlexibility::e_Bidirectional:
          {
            // assign the new sequence
            new_sse =
              util::ShPtr< assemble::SSE>
              (
                new assemble::SSE( located_sse->SubSequence( step_size, length_after_shrink), located_sse->GetType())
              );
            break;
          }
          case biol::AASequenceFlexibility::e_NTerminal:
          {
            new_sse =
              util::ShPtr< assemble::SSE>
              (
                new assemble::SSE( located_sse->SubSequence( step_size, length_after_shrink), located_sse->GetType())
              );
            break;
          }
          case biol::AASequenceFlexibility::e_CTerminal:
          {
            new_sse =
              util::ShPtr< assemble::SSE>
              (
                new assemble::SSE( located_sse->SubSequence( 0, length_after_shrink), located_sse->GetType())
              );
            break;
          }
          default:
          {
            // return empty result
            return empty_result;
          }
        }
      }

      // if extending
      else
      {
        // create a copy of the located SSE and store it in new_sse
        new_sse = util::ShPtr< assemble::SSE>( located_sse->Clone());

        // create a reference to the sequence
        const biol::AASequence &full_sequence( *PROTEIN_MODEL.GetChain( located_sse->GetChainID())->GetSequence());

        // if extending from both sides
        if( m_Side == biol::AASequenceFlexibility::e_Bidirectional || m_Side == biol::AASequenceFlexibility::e_NTerminal)
        {
          // if there are not enough residues to extend to front ( < step_size), then extend as much as you can
          const size_t extend_size( std::min< size_t>( step_size, located_sse->GetFirstAA()->GetSeqID() - 1));

          // get the residues from the full sequence
          new_sse->PrependSequence
          (
            full_sequence.SubSequence( located_sse->GetFirstAA()->GetSeqID() - extend_size - 1, extend_size)
          );
        }
        // otherwise decide on the end to extend
        if( m_Side == biol::AASequenceFlexibility::e_Bidirectional || m_Side == biol::AASequenceFlexibility::e_CTerminal)
        {
          // store number of residues in the chain
          const size_t nr_residues( full_sequence.GetSize());

          // if there are not enough residues to extend to front ( < step_size), then extend as much as you can
          const size_t extend_size( std::min( nr_residues - located_sse->GetLastAA()->GetSeqID(), step_size));

          // get the residues from the full sequence
          new_sse->AppendSequence( full_sequence.SubSequence( located_sse->GetLastAA()->GetSeqID(), extend_size));
        }

        // make sure the extending does not make this SSE overlap with any other in the protein model
        // except the original SSE we started with
        const util::SiPtrList< const assemble::SSE> overlapping_sses( PROTEIN_MODEL.GetOverlappingSSEs( *new_sse));

        // if more than one
        if( overlapping_sses.GetSize() > 1)
        {
          // warn user
          BCL_MessageVrb
          (
            "Extending causes overlap with " + util::Format()( overlapping_sses.GetSize()) + " SSEs, therefore skipping"
          );
          // return empty result
          return empty_result;
        }
      }

      // print out sse after extend/shrink
      BCL_MessageDbg( "sse after extended/shrunk: " + new_sse->GetIdentification());

      // copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // if only one end is changed, recenter the new sse
      if( m_RecenterAfterResize && m_Side != biol::AASequenceFlexibility::e_Bidirectional)
      {
        const linal::Vector3D translation_to_old_center( located_sse->GetCenter() - new_sse->GetCenter());
        new_sse->Translate( translation_to_old_center);
      }

      // insert the new SSE by replacing the old one
      if( new_sse->GetSize() > 0)
      {
        new_model->ReplaceWithOverlapping( new_sse);
      }
      else
      {
        new_model->Remove( *located_sse);
      }

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateProteinModelSSEResize::GetAlias() const
    {
      static const std::string s_alias( "MutateProteinModelSSEResize");
      return s_alias;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateProteinModelSSEResize::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Resizes SSEs in a protein model");
      serializer.AddInitializer
      (
        "locator",
        "locator for finding SSEs",
        io::Serialization::GetAgent( &m_SSELocator)
      );
      serializer.AddInitializer
      (
        "extend probability",
        "probability to extend an SSE",
        io::Serialization::GetAgent( &m_ExtendProbability)
      );
      serializer.AddInitializer
      (
        "length range",
        "range or residues to add or remove",
        io::Serialization::GetAgent( &m_LengthChangeRange)
      );
      serializer.AddInitializer
      (
        "sequence direction",
        "sequence direction of the resize",
        io::Serialization::GetAgent( &m_Side)
      );
      serializer.AddInitializer
      (
        "recenter",
        "recenter SSE after the resize",
        io::Serialization::GetAgent( &m_RecenterAfterResize)
      );
      serializer.AddInitializer
      (
        "min sse sizes",
        "minimum sizes of SSEs considered for resizing",
        io::Serialization::GetAgent( &m_MinSSESizes)
      );

      return serializer;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
