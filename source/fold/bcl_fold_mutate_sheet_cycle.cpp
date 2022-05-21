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
#include "fold/bcl_fold_mutate_sheet_cycle.h"

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
    const util::SiPtr< const util::ObjectInterface> MutateSheetCycle::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSheetCycle())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from applicable bools and ranges
    //! @param PRESERVE_SHEET_GEOMETRY boolean to whether to keep the sheet geometry intact
    //! @param PRESERVE_STRAND_ORIENTATIONS boolean to whether to keep the orientations of individiual strands intact
    //! @param ROTATE_SUBSET boolean to whether to rotate a subset
    //! @param NUMBER_ROTATIONS minimum and maximum number of rotations
    //! @param NUMBER_ROTATIONS minimum and maximum number of strands in the rotation group
    MutateSheetCycle::MutateSheetCycle
    (
      const bool PRESERVE_SHEET_GEOMETRY,
      const bool PRESERVE_STRAND_ORIENTATIONS,
      const bool ROTATE_SUBSET,
      const math::Range< size_t> &NUMBER_ROTATIONS,
      const math::Range< size_t> &SUBSET_SIZE
    ) :
      m_PreserveSheetGeometry( PRESERVE_SHEET_GEOMETRY),
      m_PreserveStrandOrientations( PRESERVE_STRAND_ORIENTATIONS),
      m_RotateSubset( ROTATE_SUBSET),
      m_NumberRotations( NUMBER_ROTATIONS),
      m_SubsetSize( SUBSET_SIZE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateSheetCycle
    MutateSheetCycle *MutateSheetCycle::Clone() const
    {
      return new MutateSheetCycle( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateSheetCycle::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes a Sheet and return a mutated Sheet
    //! @param SHEET Sheet which will be mutated
    //! @return MutateResult with the mutated Sheet
    math::MutateResult< assemble::Domain> MutateSheetCycle::operator()( const assemble::Domain &SHEET) const
    {
      // initialize an empty Sheet ShPtr to
      static util::ShPtr< assemble::Domain> s_empty_sheet;

      // make sure the passed domain has a valid topology and is of type sheet or beta-barrel
      if
      (
        !SHEET.GetTopology().IsDefined() ||
        !
        (
          SHEET.GetTopology()->GetType() == assemble::Topology::e_Sheet ||
          SHEET.GetTopology()->GetType() == assemble::Topology::e_BetaBarrel
        )
      )
      {
        // warn user and return
        BCL_MessageVrb( "The given domain is not a sheet or a barrel");
        return math::MutateResult< assemble::Domain>( s_empty_sheet, *this);
      }

      // store the number of strands
      const size_t nr_strands( SHEET.GetNumberSSEs());

      // make sure the Sheet has at least two strands
      if( nr_strands < 2)
      {
        // warn user and return
        BCL_MessageVrb( "The given sheet has less than 2 strands, therefore skipping");
        return math::MutateResult< assemble::Domain>( s_empty_sheet, *this);
      }

      // make a copy of the Sheet
      util::ShPtr< assemble::Domain> new_sheet( SHEET.Clone());

      // decide on the rotation direction
      bool rotate_direction( random::GetGlobalRandom().Boolean());

      // initialize rotation count
      const size_t nr_rotations
      (
        random::GetGlobalRandom().SizeT
        (
          math::Range< size_t>
          (
            nr_strands > ( m_NumberRotations.GetMin() + 1) ? m_NumberRotations.GetMin() : 1,
            nr_strands > ( m_NumberRotations.GetMax() + 1) ? m_NumberRotations.GetMax() : nr_strands - 1
          )
        )
      );

      // construct a vector of sses to be rotated
      util::SiPtrVector< const assemble::SSEGeometryInterface> strands_to_rotate;

      // if all of the strands will be rotated
      // or the number of strands  is already equal to/less than the given minimum subset size
      if( !m_RotateSubset || nr_strands <= m_SubsetSize.GetMin())
      {
        strands_to_rotate = SHEET.GetTopology()->GetElements();
      }
      // if a subset will be rotated
      else
      {
        // decide on the size of the subset
        const size_t subset_size
        (
          // subset size can't be larger than nr_strands
          random::GetGlobalRandom().SizeT
          (
            math::Range< size_t>
            (
              m_SubsetSize.GetMin(),
              m_SubsetSize.GetMax() >= nr_strands ? nr_strands : m_SubsetSize.GetMax()
            )
          )
        );

        // now pick a random location in the the sheet vector
        size_t index( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, nr_strands - 1)));

        // while the subset is still not full
        while( strands_to_rotate.GetSize() < subset_size)
        {
          // insert this SSE into the subset
          strands_to_rotate.PushBack( SHEET.GetTopology()->GetElements()( index));

            // increment the index
          ++index;

          // if the last is reached than set index to 0
          if( index == nr_strands)
          {
            index = 0;
          }
        }
      }

      // initialize nr_selected_strands
      const size_t nr_selected_strands( strands_to_rotate.GetSize());

      // iterate over each SSE
      for( size_t sse_index( 0); sse_index < nr_selected_strands; ++sse_index)
      {
        // calculate the index of the reference SSE, the SSE to which's position this SSE will move to
        const size_t sse_index_new
        (
          rotate_direction ?
            ( sse_index + nr_rotations) % nr_selected_strands :
            ( nr_selected_strands + sse_index - nr_rotations) % nr_selected_strands
        );

        // make a pointer of the old and reference SSE
        util::SiPtr< const assemble::SSEGeometryInterface> reference_sse( strands_to_rotate( sse_index_new));
        // make a copy of the SSE
        // this is doing dynamic cast and we have to make sure
        util::ShPtr< assemble::SSE> new_sse( strands_to_rotate( sse_index)->Clone());
        BCL_Assert( new_sse.IsDefined(), "The dynamic cast from SSEGeometryInterface to SSE failed!");

        // build up the transformation
        math::TransformationMatrix3D transformation;

        // if preserving sheet geometry
        if( m_PreserveSheetGeometry)
        {
          // first move the new SSE to origin
          transformation( math::Inverse( new_sse->GetOrientation()));

          // if the strand orientations are to be preserved
          if( m_PreserveStrandOrientations)
          {
            // calculate z twist angle
            const double z_twist
            (
              linal::Dihedral
              (
                new_sse->GetMainAxis().GetEndPoint(), new_sse->GetCenter(),
                reference_sse->GetCenter(), reference_sse->GetMainAxis().GetEndPoint()
              )
            );

            // calculate x twist angle
            const double x_twist
            (
              linal::Dihedral
              (
                new_sse->GetCenter() + new_sse->GetAxis( coord::GetAxes().e_X),
                new_sse->GetCenter(),
                reference_sse->GetCenter(),
                reference_sse->GetCenter() + reference_sse->GetAxis( coord::GetAxes().e_X)
              )
            );

            // calculate x and z orientations
            const assemble::SSEGeometryPacking::Orientation z_orientation
            (
              assemble::SSEGeometryPacking::OrientationFromTwistAngle( z_twist)
            );
            const assemble::SSEGeometryPacking::Orientation x_orientation
            (
              assemble::SSEGeometryPacking::OrientationFromTwistAngle( x_twist)
            );

            BCL_MessageDbg
            (
              "Comparing " + new_sse->GetIdentification() + " and " + reference_sse->GetIdentification()
            );
            BCL_MessageDbg
            (
              "Z-orient: " + util::Format()( math::Angle::Degree( z_twist)) + " => " + util::Format()( z_orientation)
            );
            BCL_MessageDbg
            (
              "X-orient: " + util::Format()( math::Angle::Degree( x_twist)) + " => " + util::Format()( x_orientation)
            );

            // if both z and x orientation are out of place
            if
            (
              z_orientation == assemble::SSEGeometryPacking::e_AntiParallel &&
              x_orientation == assemble::SSEGeometryPacking::e_AntiParallel
            )
            {
              BCL_MessageDbg( "Y-flip");
              // apply y flip
              transformation( coord::GetAxes().e_Y, math::g_Pi);
            }
            // if only x orientation is antiparallel
            if
            (
              z_orientation == assemble::SSEGeometryPacking::e_Parallel &&
              x_orientation == assemble::SSEGeometryPacking::e_AntiParallel
            )
            {
              BCL_MessageDbg( "Z-flip");
              // apply a z flip
              transformation( coord::GetAxes().e_Z, math::g_Pi);
            }
            // if only z orientation is antiparallel
            if
            (
              z_orientation == assemble::SSEGeometryPacking::e_AntiParallel &&
              x_orientation == assemble::SSEGeometryPacking::e_Parallel
            )
            {
              BCL_MessageDbg( "X-flip");
              // apply a x flip
              transformation( coord::GetAxes().e_X, math::g_Pi);
            }
          }

          // apply the orientation of the reference SSE
          transformation( reference_sse->GetOrientation());

        }
        // else if the sheet geometry is not to be preserved
        else
        {
          // then the transformation is just a translation
          // no need to bother preserving strand orientations, since that's done implicitly
          transformation( reference_sse->GetCenter() - new_sse->GetCenter());
        }

        // apply the transformation
        new_sse->Transform( transformation);

        // insert into the new Sheet
        new_sheet->Replace( new_sse);
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
    std::istream &MutateSheetCycle::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_PreserveSheetGeometry, ISTREAM);
      io::Serialize::Read( m_PreserveStrandOrientations, ISTREAM);
      io::Serialize::Read( m_RotateSubset, ISTREAM);
      io::Serialize::Read( m_NumberRotations, ISTREAM);
      io::Serialize::Read( m_SubsetSize, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateSheetCycle::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // read members
      io::Serialize::Write( m_PreserveSheetGeometry, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PreserveStrandOrientations, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RotateSubset, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberRotations, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SubsetSize, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
