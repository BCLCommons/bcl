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
#include "fold/bcl_fold_mutate_sheet_register_fix.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "biol/bcl_biol_atom.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    // instantiate instance
    const util::SiPtr< const util::ObjectInterface> MutateSheetRegisterFix::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSheetRegisterFix())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new MutateSheetRegisterFix
    MutateSheetRegisterFix *MutateSheetRegisterFix::Clone() const
    {
      return new MutateSheetRegisterFix( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateSheetRegisterFix::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes a Sheet and return a mutated Sheet
    //! @param SHEET Sheet which will be mutated
    //! @return MutateResult with the mutated Sheet
    math::MutateResult< assemble::Domain> MutateSheetRegisterFix::operator()( const assemble::Domain &SHEET) const
    {
      // initialize an empty Sheet ShPtr to
      static util::ShPtr< assemble::Domain> s_empty_sheet;

      // make sure the passed domain has a valid topology and is of type sheet or beta-barrel
      if
      (
        !SHEET.GetTopology().IsDefined() ||
        !(
           SHEET.GetTopology()->GetType() == assemble::Topology::e_Sheet ||
           SHEET.GetTopology()->GetType() == assemble::Topology::e_BetaBarrel
         )
      )
      {
        // warn user and return
        BCL_MessageVrb( "The given domain is not a sheet or a barrel");
        return math::MutateResult< assemble::Domain>( s_empty_sheet, *this);
      }

      // make sure the Sheet has at least two strands
      if( SHEET.GetNumberSSEs() < 2)
      {
        // warn user and return
        BCL_MessageVrb( "The given sheet has less than 2 strands, therefore skipping");
        return math::MutateResult< assemble::Domain>( s_empty_sheet, *this);
      }

      // make a copy of the Sheet
      util::ShPtr< assemble::Domain> new_sheet( SHEET.Clone());

      // get all SSEs by doing dynamic cast on the SSEGeometryInterface vector
      util::SiPtrVector< const assemble::SSE> sse_vector( SHEET.GetTopology()->GetElements());
      // check the dynamic cast
      BCL_Assert( sse_vector.IsDefined(), "The dynamic cast has failed from Topology geometries to SSEs");

      // build up the translation vector
      // first SSE stays the same
      storage::Vector< linal::Vector3D> translation_vector;
      linal::Vector3D prev_translation;
      BCL_MessageVrb( "Calculating z_translations")

      // iterate over SSEs in this sheet using the SSE vector
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
        strand_itr( sse_vector.Begin() + 1),
          strand_itr_end( sse_vector.End());
        strand_itr != strand_itr_end; ++strand_itr
      )
      {
        // create references on the ShPtrs for this strand and next strand
        const util::SiPtr< const assemble::SSE> &sp_prev_strand( *( strand_itr - 1));
        const util::SiPtr< const assemble::SSE> &sp_this_strand( *strand_itr);

        // create a reference on the SSEGeometryPacking for these two strands
        const assemble::SSEGeometryPacking &strand_pack
        (
          SHEET.GetTopology()->GetPackingForSSEGeometryPair( *sp_prev_strand, *sp_this_strand)
        );

        // find the possible hydrogen bonding residues
        const linal::Vector3D this_translation
        (
          CalculateTranslationForHydrogenBonding
          (
            strand_pack, *sp_prev_strand, *sp_this_strand
          )
        );

        BCL_MessageVrb
        (
          "sses: " + sp_prev_strand->GetIdentification() + " vs " + sp_this_strand->GetIdentification()
        );
        BCL_MessageVrb
        (
          "frags: " + strand_pack.GetFirstSSEGeometry()->GetIdentification() +
          " vs " + strand_pack.GetSecondSSEGeometry()->GetIdentification()
        );
        BCL_MessageVrb
        (
          "translation: " + util::Format()( this_translation)
        );

        // update the previous translation
        prev_translation += this_translation;

        BCL_MessageVrb
        (
          "prev translation after update: " + util::Format()( prev_translation)
        );

        // update the vector
        translation_vector.PushBack( prev_translation);
      }

      // create iterator on the z_translation_vector
      storage::Vector< linal::Vector3D>::const_iterator translation_itr( translation_vector.Begin());
      const storage::Vector< linal::Vector3D>::const_iterator translation_itr_end( translation_vector.End());

      BCL_MessageVrb( "Appyling z_translations");

      // now iterate again over the Sheet to apply the translations
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
        strand_itr( sse_vector.Begin() + 1),
          strand_itr_end( sse_vector.End());
        strand_itr != strand_itr_end && translation_itr != translation_itr_end;
        ++strand_itr, ++translation_itr
      )
      {
        // create references on the ShPtrs for this strand and next strand
        const util::SiPtr< const assemble::SSE> &sp_prev_strand( *( strand_itr - 1));
        const util::SiPtr< const assemble::SSE> &sp_this_strand( *strand_itr);

        // create a reference on the SSEGeometryPacking for these two strands
        const assemble::SSEGeometryPacking &strand_pack
        (
          SHEET.GetTopology()->GetPackingForSSEGeometryPair( *sp_prev_strand, *sp_this_strand)
        );

        BCL_MessageVrb
        (
          "sses: " + sp_prev_strand->GetIdentification() + " vs " + sp_this_strand->GetIdentification()
        );
        BCL_MessageVrb
        (
          "frags: " + strand_pack.GetFirstSSEGeometry()->GetIdentification() +
          " vs " + strand_pack.GetSecondSSEGeometry()->GetIdentification()
        );
        BCL_MessageVrb( "translation: " + util::Format()( *translation_itr));

        // make a copy of this strand
        util::ShPtr< assemble::SSE> new_strand( sp_this_strand->Clone());

        // apply the translation
        new_strand->Translate( *translation_itr);

        // do replace
        new_sheet->Replace( new_strand);
      }

      // end
      return math::MutateResult< assemble::Domain>( new_sheet, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateSheetRegisterFix::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateSheetRegisterFix::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief finds possible hydrogen bonding pair residues from a StrandPacking objects
    //! @param STRAND_PACK SSEGeometryPacking between two strands
    //! @param STRAND_A first strand of interest
    //! @param STRAND_B second strand of interest
    //! @return possible hydrogen bonding pair residues from a StrandPacking objects
    linal::Vector3D MutateSheetRegisterFix::CalculateTranslationForHydrogenBonding
    (
      const assemble::SSEGeometryPacking &STRAND_PACK,
      const assemble::SSE &STRAND_A,
      const assemble::SSE &STRAND_B
    )
    {
      // create references on the first and second sequence and sse pack
//      const SSE &sse_a( *STRAND_PACK.First().First());
//      const SSE &sse_b( *STRAND_PACK.First().Second());
      const assemble::SSE &sse_a( STRAND_A);
      const assemble::SSE &sse_b( STRAND_B);

      // create boolean whether the first sequence was hydrogen donator
      storage::VectorND< 2, util::SiPtr< const biol::AABase> > best_pair;
      double best_distance( 20.0);

      // first locate the closest amino acids
      // iterate over both sequences
      for
      (
        biol::AASequence::const_iterator aa_itr_a( sse_a.Begin()), aa_itr_a_end( sse_a.End());
        aa_itr_a != aa_itr_a_end; ++aa_itr_a
      )
      {
        for
        (
          biol::AASequence::const_iterator aa_itr_b( sse_b.Begin()), aa_itr_b_end( sse_b.End());
          aa_itr_b != aa_itr_b_end; ++aa_itr_b
        )
        {
          // calculate distance between N of aa_itr_a, and O of aa_itr_b
          const double this_distance
          (
            linal::Distance( ( *aa_itr_a)->GetCA().GetCoordinates(), ( *aa_itr_b)->GetCA().GetCoordinates())
          );

          // if the distance is less than the best distance
          if( this_distance < best_distance)
          {
            // find if CA-CB vector face the same direction
            // calculate dihedral
            const double dihedral
            (
              linal::Dihedral
              (
                ( *aa_itr_a)->GetFirstSidechainAtom().GetCoordinates(),
                ( *aa_itr_a)->GetCA().GetCoordinates(),
                ( *aa_itr_b)->GetCA().GetCoordinates(),
                ( *aa_itr_b)->GetFirstSidechainAtom().GetCoordinates()
              )
            );

            // if they are facing the same way
            if( dihedral < math::g_Pi / 2)
            {
              // update the best distance and pair information
              best_distance = this_distance;
              best_pair.First() = **aa_itr_a;
              best_pair.Second() = **aa_itr_b;
            }
          }
        }
      }

      BCL_MessageVrb
      (
        "CAs to adjust " + best_pair.First()->GetIdentification() + " and " + best_pair.Second()->GetIdentification()
        + " at distance " + util::Format()( best_distance)
      );

      // find the expected CA coordinate ( this will be off the z axis of sse_b)
      const linal::Vector3D ca_expected
      (
        best_pair.First()->GetCA().GetCoordinates() + STRAND_PACK.GetShortestConnection().GetDirection()
      );

      BCL_MessageVrb( "expect CA position: " + util::Format()( ca_expected));
      // now find the translation required from current_ca to ca_expected
      const linal::Vector3D translation
      (
        ( ca_expected - best_pair.Second()->GetCA().GetCoordinates())
      );

      // now find the component of this translation along the z-axis
      const double z_translation
      (
        translation * sse_b.GetAxis( coord::GetAxes().e_Z)
      );

      // return value
      return z_translation * sse_b.GetAxis( coord::GetAxes().e_Z);
    }

  } // namespace fold
} // namespace bcl
