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
#include "fold/bcl_fold_placement_strand_next_to_sheet.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PlacementStrandNextToSheet::s_Instance
    (
      util::Enumerated< PlacementInterface< assemble::SSE, assemble::ProteinModel> >::AddInstance( new PlacementStrandNextToSheet())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a flip probability and a distance deviation
    //! @param FLIP_PROBABILITY probability of applying a flip the strand around a random axis
    PlacementStrandNextToSheet::PlacementStrandNextToSheet( const double FLIP_PROBABILITY) :
      m_FlipProbability( FLIP_PROBABILITY)
    {
      BCL_Assert
      (
        m_FlipProbability >= 0.0 && m_FlipProbability <= 1.0,
        "The flip probability should be in the range [0..1] not " + util::Format()( FLIP_PROBABILITY)
      );
    }

    //! @brief Clone function
    //! @return pointer to new PlacementStrandNextToSheet
    PlacementStrandNextToSheet *PlacementStrandNextToSheet::Clone() const
    {
      return new PlacementStrandNextToSheet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &PlacementStrandNextToSheet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &PlacementStrandNextToSheet::GetAlias() const
    {
      static const std::string s_alias( "PlacementStrandNextToSheet");
      return s_alias;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PlacementStrandNextToSheet::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Places a strand next to a sheet.");
      serializer.AddInitializer
      (
        "flip_probability",
        "probability of flipping a strand around a random axis",
        io::Serialization::GetAgent( &m_FlipProbability)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief find the placement for the given strand SELECTED_STRAND with respect to a sheet in PROTEIN_MODEL
    //! @param SELECTED_STRAND SiPtr to strand to be placed
    //! @param PROTEIN_MODEL to which the SELECTED_STRAND is going to be added
    //! @return the transformationmatrix3d to place the SELECTED_STRAND in the PROTEIN_MODEL
    storage::Pair< math::TransformationMatrix3D, bool> PlacementStrandNextToSheet::Place
    (
      const assemble::SSE &SELECTED_STRAND,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize a static transformation and boolean pair to return in case of failures
      static const storage::Pair< math::TransformationMatrix3D, bool>
        s_failure_return( math::TransformationMatrix3D(), false);

      // make sure the provided sse is strand
      BCL_Assert
      (
        SELECTED_STRAND.GetType() == biol::GetSSTypes().STRAND,
        "The provided SSE should be a strand not of type " + SELECTED_STRAND.GetType().GetName()
      );

      // collect the sheets from the given ProteinModel
      util::ShPtrVector< assemble::Domain> collected_sheets( assemble::CollectorSheet().Collect( PROTEIN_MODEL));

      // initialize variable to hold the sheet of interest
      util::ShPtr< assemble::Domain> this_sheet;

      // if there are no sheets found
      if( collected_sheets.IsEmpty())
      {
        // warn user and return failure
        BCL_MessageVrb( "No sheets found in the model therefore skipping this move");
        return s_failure_return;
      }
      // if there is only one sheet than pick it
      else if( collected_sheets.GetSize() == 1)
      {
        this_sheet = collected_sheets.FirstElement();
      }
      // otherwise pick one randomly
      else
      {
        this_sheet =
          collected_sheets
          (
            random::GetGlobalRandom().SizeT
            (
              math::Range< size_t>( 0, collected_sheets.GetSize() - 1)
            )
          );
      }

      // make sure the sheet is not empty
      BCL_Assert( this_sheet->GetNumberSSEs() > 0, "The selected sheet is empty this should not have happened!!");

      // issue message
      BCL_MessageVrb( "Following sheet was selected:\n " + this_sheet->GetIdentification());

      // calculate the placement and return
      return Place( SELECTED_STRAND, *this_sheet);
    }

    //! @brief find the placement for the given strand SELECTED_STRAND with respect to a sheet in SHEET
    //! @param SELECTED_STRAND SiPtr to strand to be placed
    //! @param SHEET to which the SELECTED_STRAND is going to be added
    //! @return the transformationmatrix3d to place the SELECTED_STRAND in the SHEET
    storage::Pair< math::TransformationMatrix3D, bool> PlacementStrandNextToSheet::Place
    (
      const assemble::SSE &SELECTED_STRAND,
      const assemble::Domain &SHEET
    ) const
    {
      // initialize a static transformation and boolean pair to return in case of failures
      static const storage::Pair< math::TransformationMatrix3D, bool>
        s_failure_return( math::TransformationMatrix3D(), false);

      // make sure the passed domain has a valid topology and is of type sheet or beta-barrel
      if
      (
        !SHEET.GetTopology().IsDefined() || SHEET.GetTopology()->GetType() != assemble::Topology::e_Sheet
      )
      {
        // warn user and return
        BCL_MessageVrb( "The given domain is not a sheet");
        return s_failure_return;
      }

      // initialize the transformation matrix that will be applied to this sse
      math::TransformationMatrix3D transformation;

      // initialize the distance using the randomly determined from the given ranges
      const double distance( random::GetGlobalRandom().Double( contact::GetTypes().STRAND_STRAND->GetPreferredDistanceRange()));

      // if the sheet has only one strand
      if( SHEET.GetNumberSSEs() == 1)
      {
        BCL_MessageVrb( "only one strand found in the sheet");
        // make a reference on the only strand in the sheet
        const assemble::SSE &this_strand( *SHEET.GetSSEs().FirstElement());

        // get a random sign to determine direction
        const int direction( random::GetGlobalRandom().Sign());

        // get a random twist angle between allowed ranges and multiply with translation sign
        // since we do the rotation always around the y axis
        const double twist_angle
        (
          random::GetGlobalRandom().Double( contact::GetTypes().STRAND_STRAND->GetTiltAngleRange()) * direction
        );

        // initialize the translation vector that translate along Y axis
        linal::Vector3D this_translation( 0, direction * distance, 0);

        // rotate the translation vector with neighbor sses rotation
        this_translation.Rotate( this_strand.GetRotation());

        BCL_MessageVrb( "The twist angle is " + util::Format()( math::Angle::Degree( twist_angle)));
        // apply the twist rotation
        transformation( math::RotationMatrix3D( coord::GetAxes().e_Y, twist_angle));

        // apply the orientation of the strand
        transformation( this_strand.GetOrientation());

        // apply the calculated translation to transformation
        transformation( this_translation);
      }
      // if more than one strand
      else
      {
        // decide whether to add before or after
        const bool add_after( random::GetGlobalRandom().Boolean());

        BCL_MessageVrb( "multiple strands in the sheet, add after: " + util::Format()( add_after));

        // initialize variable to store the edge strand
        const util::SiPtr< const assemble::SSEGeometryInterface> sp_edge_strand
        (
          add_after ?
            SHEET.GetTopology()->GetElements().LastElement() :
            SHEET.GetTopology()->GetElements().FirstElement()
        );
        // initialize variables to store the edge_strand and reference_strand
        const util::SiPtr< const assemble::SSEGeometryInterface> sp_reference_strand
        (
          add_after ?
            SHEET.GetTopology()->GetElements()( SHEET.GetNumberSSEs() - 2) :
            SHEET.GetTopology()->GetElements()( 1)
        );

        // get the StrandPacking information
        const assemble::SSEGeometryPacking &this_packing
        (
          SHEET.GetTopology()->GetPackingForSSEGeometryPair( *sp_edge_strand, *sp_reference_strand)
        );

        // make sure the orientation was found
        BCL_Assert
        (
          this_packing.GetOrientation() != assemble::SSEGeometryPacking::e_UndefinedOrientation,
          "Orientation was undefined for these two strands!"
        );

        BCL_MessageVrb
        (
          "edge_strand: " + sp_edge_strand->GetIdentification() +
          " reference_strand: " + sp_reference_strand->GetIdentification()
        );
        BCL_MessageVrb
        (
          "edge_frag: " + this_packing.GetFirstSSEGeometry()->GetIdentification() +
          " reference_frag: " + this_packing.GetSecondSSEGeometry()->GetIdentification()
        );

        // get the translation vector from reference to edge
        linal::Vector3D translation_vector
        (
          this_packing.GetFirstSSEGeometry()->GetCenter() - this_packing.GetSecondSSEGeometry()->GetCenter()
        );

        // make a copy of the transformation of selected_strand and edge
        math::TransformationMatrix3D selected_trans( SELECTED_STRAND.GetOrientation());
        math::TransformationMatrix3D edge_trans( sp_edge_strand->GetOrientation());

        // calculate the twist angle from edge to reference
        double twist_angle
        (
          linal::Dihedral
          (
            sp_edge_strand->GetCenter() + sp_edge_strand->GetAxis( coord::GetAxes().e_Z),
            sp_edge_strand->GetCenter(),
            sp_reference_strand->GetCenter(),
            sp_reference_strand->GetCenter() + sp_reference_strand->GetAxis( coord::GetAxes().e_Z)
          )
        );
        // calculate the line connecting two origins
        linal::Vector3D next_spot( sp_edge_strand->GetCenter() + translation_vector);
        const double translation_x_proj_angle
        (
          linal::ProjAngle( translation_vector, sp_edge_strand->GetAxis( coord::GetAxes().e_Y))
        );

        BCL_MessageVrb
        (
          "The translation-x proj angle is: " + util::Format()( math::Angle::Degree( translation_x_proj_angle)) + " degrees "
          + util::Format()( translation_x_proj_angle) + " radians"
        );
        // if projection angle of x axis is smaller than the 90 degrees with respect to translation
        // meaning it was anti-parallel
        if( translation_x_proj_angle <= ( math::g_Pi / 2.0))
        {
          BCL_MessageVrb
          (
            "The translation-x proj angle is less than 90 degrees"
          );
          // revert the twist angle
          twist_angle *= -1.0;
        }

        // apply the twist angle rotation
        transformation( math::RotationMatrix3D( coord::GetAxes().e_Y, twist_angle));

        // apply the transformation of the edge strand
        transformation( sp_edge_strand->GetOrientation());

        // apply the translation
        transformation( translation_vector);
      }

      // if flip was selected
      if( random::GetGlobalRandom().Double() < m_FlipProbability)
      {
        // decide on around which axis to flip around ( X or Y or Z)
        const coord::Axis flip_axis( random::GetGlobalRandom().Random< size_t>( 2));

        // issue message
        BCL_MessageVrb( "Flip will be applied along axis " + flip_axis.GetName());

        // make a copy of the transformation
        math::TransformationMatrix3D transformation_copy( transformation);

        transformation( math::Inverse( transformation_copy));
        // apply a flip along the x or y axis
        transformation( flip_axis, math::g_Pi);
        transformation( transformation_copy);
      }

      // initialize the transformation matrix that will be applied to this sse
      return storage::Pair< math::TransformationMatrix3D, bool>( transformation, true);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
