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
#include "fold/bcl_fold_mutate_protein_model_sse_pair_align_and_pull.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "score/bcl_score_aa_pair_hi_res_clash.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSEPairAlignAndPull::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelSSEPairAlignAndPull())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSEPairAlignAndPull::MutateProteinModelSSEPairAlignAndPull() :
      m_Collector(),
      m_Pull( false),
      m_Scheme( GetStaticClassName< MutateProteinModelSSEPairAlignAndPull>())
    {
    }

    //! @brief construct from Collector, max translation and rotation
    //! @param SP_COLLECTOR ShPtr to collector of SSE pairs
    //! @param MAX_TRANSLATION maximum translation allowed
    //! @param MAX_ROTATION maximum rotation allowed
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEPairAlignAndPull::MutateProteinModelSSEPairAlignAndPull
    (
      const util::ShPtr
      <
        find::CollectorInterface
        <
          storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface
        >
      > &SP_COLLECTOR,
      const bool &PULL,
      const std::string &SCHEME
    ) :
      m_Collector( SP_COLLECTOR),
      m_Pull( PULL),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSEPairAlignAndPull
    MutateProteinModelSSEPairAlignAndPull *MutateProteinModelSSEPairAlignAndPull::Clone() const
    {
      return new MutateProteinModelSSEPairAlignAndPull( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelSSEPairAlignAndPull::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSEPairAlignAndPull::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // if no pair is available
      if( PROTEIN_MODEL.GetNumberSSEs() < size_t( 2))
      {
        // warn user
        BCL_MessageVrb( "unable to find a pair of sses");

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      auto sses( PROTEIN_MODEL.GetSSEs());
      size_t sse_id_a( random::GetGlobalRandom().Random( sses.GetSize() - 1));
      size_t sse_id_b( random::GetGlobalRandom().Random( sses.GetSize() - 2));
      if( sse_id_b >= sse_id_a)
      {
        ++sse_id_b;
      }
      static const double s_max_align_dist_wo_pull( 20.0);
      static const size_t s_max_tries( 100);
      for
      (
        size_t tries( 0);
        !m_Pull
        && linal::Distance( sses( sse_id_a)->GetCenter(), sses( sse_id_b)->GetCenter()) > s_max_align_dist_wo_pull
        && tries < s_max_tries;
        ++tries
      )
      {
        sse_id_a = random::GetGlobalRandom().Random( sses.GetSize() - 1);
        sse_id_b = random::GetGlobalRandom().Random( sses.GetSize() - 2);
        if( sse_id_b >= sse_id_a)
        {
          ++sse_id_b;
        }
      }

      // calculate the sse_packing
      const assemble::SSEGeometryPacking sse_pack( *sses( sse_id_a), *sses( sse_id_b));

      // pick either sse
      util::ShPtr< assemble::SSE> mutated_sse
      (
        sses( sse_id_a)->Clone()
      );

      BCL_MessageVrb( "mutated sse " + mutated_sse->GetIdentification() + " relative to " + sses( sse_id_b)->GetIdentification());

      // generate TransformationMatrix3D and place sse in origin
      const linal::Vector3D center_a( mutated_sse->GetCenter());
      mutated_sse->Translate( -center_a);
      coord::LineSegment3D main_axis_recentered_b
      (
        ( sses( sse_id_b)->GetMainAxis().GetStartPoint() - sses( sse_id_b)->GetMainAxis().GetEndPoint()) / 2.0,
        ( sses( sse_id_b)->GetMainAxis().GetEndPoint() - sses( sse_id_b)->GetMainAxis().GetStartPoint()) / 2.0
      );
      if( random::GetGlobalRandom().Boolean())
      {
        main_axis_recentered_b = coord::LineSegment3D( main_axis_recentered_b.GetStartPoint(), main_axis_recentered_b.GetEndPoint());
      }
      math::TransformationMatrix3D transform( mutated_sse->GetMainAxis(), main_axis_recentered_b);
      mutated_sse->Transform( transform);
      mutated_sse->Translate( center_a);

      linal::Vector3D translation_dir( sses( sse_id_b)->GetCenter() - center_a);
      static score::AAPairHiResClash s_aaclash;
      const double true_dist( translation_dir.Norm());
      translation_dir.Normalize();
      transform.SetTranslation( linal::Vector3D());

      mutated_sse->Transform( transform);
      if( m_Pull && s_aaclash( *mutated_sse, *sses( sse_id_b)) < 0.5)
      {
        double max_remaining_translation( std::max( true_dist - 5.0, 0.0));
        while( max_remaining_translation > 0.1)
        {
          double amnt_to_translate( std::max( 0.1, max_remaining_translation * 0.5));
          mutated_sse->Translate( translation_dir * amnt_to_translate);
          if( s_aaclash( *mutated_sse, *sses( sse_id_b)) > 0.5)
          {
            mutated_sse->Translate( translation_dir * -amnt_to_translate);
            max_remaining_translation = amnt_to_translate;
          }
          else
          {
            max_remaining_translation -= amnt_to_translate;
          }
        }
      }

      // make copy of proteinmodel
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // replace mutated sse
      new_model->Replace( mutated_sse);

      // return
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelSSEPairAlignAndPull::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Collector, ISTREAM);
      io::Serialize::Read( m_Pull, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateProteinModelSSEPairAlignAndPull::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Collector, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Pull, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
