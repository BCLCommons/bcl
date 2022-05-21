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
#include "fold/bcl_fold_mutate_sheet_register_shift.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateSheetRegisterShift::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSheetRegisterShift())
    );

    //! @brief static function to return default shift probabilities
    //! @return default shift probabilities
    const storage::VectorND< 2, double> &MutateSheetRegisterShift::GetDefaultShiftProbabilities()
    {
      // initialize static const shift probabilities vector
      static const storage::VectorND< 2, double> s_default_probabilities( 0.5, 0.5);

      // end
      return s_default_probabilities;
    }

    //! @brief static function to return CA z-translation for a single residue
    //! @return CA z-translation for a single residue
    double MutateSheetRegisterShift::GetCAShiftAlongZAxis()
    {
      // have a static double
      static const double s_ca_shift( 3.425);

      // end
      return s_ca_shift;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a shift probabilities vector
    //! @param SHIFT_PROBABILITIES vector that contains probabilities from regular and flip shift
    MutateSheetRegisterShift::MutateSheetRegisterShift
    (
      const storage::VectorND< 2, double> &SHIFT_PROBABILITIES
    ) :
      m_ShiftProbabilities( SHIFT_PROBABILITIES)
    {
      // make sure the probabilities are both positive
      BCL_Assert
      (
        m_ShiftProbabilities.First() >= 0.0 && m_ShiftProbabilities.Second() >= 0.0,
        "The shift probabilities should be both >= 0.0 not " + util::Format()( m_ShiftProbabilities)
      );

      // normalize the probabilities
      const double sum( m_ShiftProbabilities.First() + m_ShiftProbabilities.Second());
      m_ShiftProbabilities.First() /= sum;
      m_ShiftProbabilities.Second() /= sum;
    }

    //! @brief Clone function
    //! @return pointer to new MutateSheetRegisterShift
    MutateSheetRegisterShift *MutateSheetRegisterShift::Clone() const
    {
      return new MutateSheetRegisterShift( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateSheetRegisterShift::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes a Sheet and return a mutated Sheet
    //! @param SHEET Sheet which will be mutated
    //! @return MutateResult with the mutated Sheet
    math::MutateResult< assemble::Domain> MutateSheetRegisterShift::operator()( const assemble::Domain &SHEET) const
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

      // get a random iterator
      util::SiPtrVector< const assemble::SSE>::const_iterator prev_strand_itr
      (
        random::GetGlobalRandom().Iterator
        (
          sse_vector.Begin(),
          sse_vector.End() - 1,
          sse_vector.GetSize() - 1
        )
      );

      // create an iterator to store the reference sse
      util::SiPtrVector< const assemble::SSE>::const_iterator this_strand_itr( prev_strand_itr + 1);

      // create references on the ShPtrs for this strand and next strand
      const util::SiPtr< const assemble::SSE> &sp_prev_strand( *prev_strand_itr);
      const util::SiPtr< const assemble::SSE> &sp_this_strand( *this_strand_itr);

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

      // determine whether to shift or flip-twist
      bool flip_shift( random::GetGlobalRandom().Double() >= m_ShiftProbabilities.First());

      // determine random direction for the shift
      int direction( random::GetGlobalRandom().Sign());

      // initialize translation
      linal::Vector3D translation;

      // make a copy of this SSE
      util::ShPtr< assemble::SSE> new_strand( sp_this_strand->Clone());

      if( flip_shift)
      {
        BCL_MessageVrb( "applying flip shift");
        // set the translation along the Z axis in random direction of length from first CA to the second CA
        translation =
          strand_pack.GetSecondSSEGeometry()->GetAxis( coord::GetAxes().e_Z) * direction * GetCAShiftAlongZAxis();

        // build up the z flip transformation
        math::TransformationMatrix3D transformation( math::Inverse( sp_this_strand->GetOrientation()));
        transformation( coord::GetAxes().e_Z, math::g_Pi);
        transformation( sp_this_strand->GetOrientation());

        // apply the transformation
        new_strand->Transform( transformation);
      }
      // if no flip
      else
      {
        BCL_MessageVrb( "applying regular shift");
        // set the translation along the Z axis in random direction of length from first CA to third CA
        translation =
          strand_pack.GetSecondSSEGeometry()->GetAxis( coord::GetAxes().e_Z) * direction * GetCAShiftAlongZAxis() * 2;
      }

      // apply the translation
      new_strand->Translate( translation);

      BCL_MessageVrb( "translation: " + util::Format()( translation));

      // replace into new sheet
      new_sheet->Replace( new_strand);

      // now iterate over the strands after
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          strand_itr( this_strand_itr + 1), strand_itr_end( sse_vector.End());
        strand_itr != strand_itr_end;
        ++strand_itr
      )
      {
        // make a hardcopy of this SSE
        util::ShPtr< assemble::SSE> copy_strand( ( *strand_itr)->Clone());

        // if the flip_shift was selected
        if( flip_shift)
        {
          // build up the z flip transformation
          math::TransformationMatrix3D transformation( math::Inverse( copy_strand->GetOrientation()));
          transformation( coord::GetAxes().e_Z, math::g_Pi);
          transformation( copy_strand->GetOrientation());

          // apply the transformation
          copy_strand->Transform( transformation);
        }

        // apply the translation
        copy_strand->Translate( translation);

        // replace into new sheet
        new_sheet->Replace( copy_strand);
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
    std::istream &MutateSheetRegisterShift::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ShiftProbabilities, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateSheetRegisterShift::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ShiftProbabilities, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
