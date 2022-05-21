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
#include "fold/bcl_fold_placement_sse_short_loop.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "contact/bcl_contact_types.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! minimum loop distance
    const double PlacementSSEShortLoop::s_MinLoopDistance( 2.0);

    //! maxiumum loop distance
    const double PlacementSSEShortLoop::s_MaxLoopDistance( 15.0);

    //! maximum loop distance per residue
    const double PlacementSSEShortLoop::s_MaxDistancePerResidue( 3.04964);

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PlacementSSEShortLoop::s_Instance
    (
      util::Enumerated< PlacementInterface< assemble::SSE, assemble::ProteinModel> >::AddInstance( new PlacementSSEShortLoop())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PlacementSSEShortLoop::PlacementSSEShortLoop() :
      m_MaxShortLoopLength( 7),
      m_AddToTopProbability( 0.5),
      m_MaxHingeAngleAddingToTop( math::g_Pi / 3.0)
    {
    }

    //! @brief constructor from a maximum short loop length
    //! @param MAX_SHORT_LOOP_LENGTH maximum number of residues between SSEs to be classified as short loop
    //! @param ADD_TO_TOP_PROBABILITY probability that determines how frequently the SSE will be added to the top of the SSE
    //! @param MAX_HINGE_ANGLE_ADDING_TO_TOP max hinge angle to be used when adding to top
    PlacementSSEShortLoop::PlacementSSEShortLoop
    (
      const size_t MAX_SHORT_LOOP_LENGTH,
      const double ADD_TO_TOP_PROBABILITY,
      const double MAX_HINGE_ANGLE_ADDING_TO_TOP
    ) :
      m_MaxShortLoopLength( MAX_SHORT_LOOP_LENGTH),
      m_AddToTopProbability( ADD_TO_TOP_PROBABILITY),
      m_MaxHingeAngleAddingToTop( MAX_HINGE_ANGLE_ADDING_TO_TOP)
    {

    }

    //! @brief Clone function
    //! @return pointer to new PlacementSSEShortLoop
    PlacementSSEShortLoop *PlacementSSEShortLoop::Clone() const
    {
      return new PlacementSSEShortLoop( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PlacementSSEShortLoop::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate placement for given SSE at a random orientation wrt to a located SSE in PROTEIN_MODEL
    //! @param SELECTED_SSE SiPtr to SSE to be placed, it needs to be in the origin
    //! @param PROTEIN_MODEL to which the SSE is going to be added
    //! @return the transformationmatrix3d to place the SSE in the PROTEIN_MODEL
    storage::Pair< math::TransformationMatrix3D, bool> PlacementSSEShortLoop::Place
    (
      const assemble::SSE &SELECTED_SSE,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // store eligible sses that have short loops to SELECTED_SSE
      util::SiPtrList< const assemble::SSE> eligible_sses
      (
        PROTEIN_MODEL.GetSSEsWithShortLoops( SELECTED_SSE, m_MaxShortLoopLength)
      );

      // initialize neighbor_sse
      util::SiPtr< const assemble::SSE> neighbor_sse;

      // if no eligible sses are found
      if( eligible_sses.IsEmpty())
      {
        // message
        BCL_MessageVrb( "No sses with short loops were found!");

        // return false placement
        return storage::Pair< math::TransformationMatrix3D, bool>( math::TransformationMatrix3D(), false);
      }

      // else if there is only one eligible sse
      if( eligible_sses.GetSize() == 1)
      {
        // return that sses transformation
        neighbor_sse = eligible_sses.FirstElement();
      }
      // if more than one pick one
      else
      {
        neighbor_sse = assemble::PickSSERandom().Pick( eligible_sses);
      }

      // if the neighbor sse is empty
      if( !neighbor_sse.IsDefined())
      {
        // give a warning and return
        BCL_MessageCrt( "SSE picked with short loop is undefined!!!");
        return storage::Pair< math::TransformationMatrix3D, bool>( math::TransformationMatrix3D(), false);
      }

      // calculate the placement and return it
      return Place( SELECTED_SSE, *neighbor_sse);
    }

    //! @brief generate placement for given SSE at a random orientation wrt to NEIGHBOR_SSE
    //! @param SELECTED_SSE SiPtr to SSE to be placed
    //! @param NEIGHBOR_SSE neighbor SSE
    //! @return the transformationmatrix3d to place the SSE next to NEIGHBOR_SSE
    storage::Pair< math::TransformationMatrix3D, bool> PlacementSSEShortLoop::Place
    (
      const assemble::SSE &SELECTED_SSE,
      const assemble::SSE &NEIGHBOR_SSE
    ) const
    {
      // if any of the given SSEs are loops then return undefined transformation
      if( !SELECTED_SSE.GetType()->IsStructured() || !NEIGHBOR_SSE.GetType()->IsStructured())
      {
        return storage::Pair< math::TransformationMatrix3D, bool>( math::TransformationMatrix3D(), false);
      }

      // create the contact type
      contact::Type contact_type( contact::GetTypes().TypeFromSSTypes( SELECTED_SSE, NEIGHBOR_SSE));

      // determine if selected_sse comes before or after neighbor sse
      const bool selected_sse_comes_first( SELECTED_SSE.GetFirstAA()->GetSeqID() < NEIGHBOR_SSE.GetFirstAA()->GetSeqID());

      // calculate the sequence distance and loop distance
      const size_t seq_dist( biol::CalculateSequenceDistance( SELECTED_SSE, NEIGHBOR_SSE));
      const double ideal_loop_distance( s_MaxDistancePerResidue * seq_dist);
      const double max_loop_distance( std::min( s_MaxLoopDistance, ideal_loop_distance));

      // decide whether to to regular placement or add at the end
      const bool add_at_top( random::GetGlobalRandom().Double() < m_AddToTopProbability);

      // if adding at the top
      if( add_at_top)
      {
        // calculate maximum allowed loop distance
        const double loop_distance
        (
          seq_dist == 0 ?
          s_MinLoopDistance : // if broken SSEs use min loop distance, otherwise use random within range
          random::GetGlobalRandom().Double( math::Range< double>( s_MinLoopDistance, max_loop_distance))
        );

        // calculate z translation
        linal::Vector3D z_translation
        (
          0.0, 0.0,
          double( math::ConvertBooleanToSign( !selected_sse_comes_first)) *
          ( ( NEIGHBOR_SSE.GetLength() / 2.0) + ( SELECTED_SSE.GetLength() / 2.0) + loop_distance)
        );

        // decide on which direction the SSE should be facing
        const double arc_angle( random::GetGlobalRandom().Random< double>( 0, 2.0 * math::g_Pi));

        // determine rotation axis
        linal::Vector3D rotation_axis( cos( arc_angle), sin( arc_angle), 0.0);

        // decide on a hinge angle
        // 0 is complete z translation only, while -90 or 90 means perpendicular
        const double hinge_angle
        (
          random::GetGlobalRandom().Double
          (
            math::Range< double>( -m_MaxHingeAngleAddingToTop, m_MaxHingeAngleAddingToTop)
          )
        );

        // form transformation
        math::TransformationMatrix3D transformation;

        // apply a random rotation around z axis
        const double z_rotation_angle( random::GetGlobalRandom().Random< double>( 0.0, 2.0 * math::g_Pi));
        transformation( coord::GetAxes().e_Z, z_rotation_angle);

        // find the translation to move the end of the SSE which will connect the loop to origin
        // if selected SSE comes first it's the C terminal edge, so we have to move down
        // otherwise it's the N terminal edge, so we have to move up
        linal::Vector3D translate_loop_end_to_center
        (
          0.0, 0.0, SELECTED_SSE.GetLength() / 2.0 * double( math::ConvertBooleanToSign( !selected_sse_comes_first))
        );

        // move the loop end to origin
        transformation( translate_loop_end_to_center);

        // apply the tilt angle
        transformation( math::RotationMatrix3D( rotation_axis, hinge_angle));

        // move it back up
        transformation( -translate_loop_end_to_center);

        // apply translation
        transformation( z_translation);

        // apply the neighbor sse's orientation
        transformation( NEIGHBOR_SSE.GetOrientation());

        // return
        return storage::Pair< math::TransformationMatrix3D, bool>( transformation, true);

      }
      // adding next to it
      else
      {
        // decide on a flip axis
        const coord::Axis flip_axis( random::GetGlobalRandom().Boolean() ? coord::GetAxes().e_X : coord::GetAxes().e_Y);

        // initialize translation direction
        linal::Vector3D translation_direction;

        // initialize the transformation matrix that will be applied to this sse
        math::TransformationMatrix3D transformation;

        // if helix helix
        if( contact_type == contact::GetTypes().HELIX_HELIX)
        {
          // determine a random angle between -90 and +90
          const double angle( random::GetGlobalRandom().Random< double>( 0.0, 2.0 * math::g_Pi));

          // use this angle to determine a random unit vector on xy plane and scale the translation
          translation_direction = linal::Vector3D( cos( angle), sin( angle), 0.0);
        }
        // if the selected_sse is a helix and the neighbor sse is strand
        else if( contact_type == contact::GetTypes().HELIX_SHEET)
        {
          // then we can only use x axis
          translation_direction = linal::Vector3D( random::GetGlobalRandom().Sign(), 0.0, 0.0);
        }
        // if the selected SSE is a strand and the neighbor sse is helix
        else if( contact_type == contact::GetTypes().SHEET_HELIX)
        {
          // determine a random angle between -90 and +90
          const double angle( random::GetGlobalRandom().Random< double>( 0.0, 2.0 * math::g_Pi));

          // use this angle to determine a random unit vector on xy plane
          translation_direction = linal::Vector3D( cos( angle), sin( angle), 0.0);

          // we also need a rotation around Z to make sure the strand side chains are facing the helix
          transformation( coord::GetAxes().e_Z, -angle);

        }
        else if( contact_type == contact::GetTypes().STRAND_STRAND)
        {
          // determine randomly to use x (SHEET_SHEET) or y axis( STRAND_STRAND)
          const coord::Axis trans_axis( random::GetGlobalRandom().Boolean() ? coord::GetAxes().e_X : coord::GetAxes().e_Y);

          // set the translation direction
          translation_direction( trans_axis) = random::GetGlobalRandom().Sign();

          // if sheet sheet was selected update the distance and tilt angle
          if( trans_axis == coord::GetAxes().e_X)
          {
            // update contact_type
            contact_type = contact::GetTypes().SHEET_SHEET;
          }
        }
        else
        {
          // if no valid contact type found then return undefined transformation
          BCL_MessageCrt( "undefined contact type!" + contact_type.GetName());
          return storage::Pair< math::TransformationMatrix3D, bool>( math::TransformationMatrix3D(), false);
        }

        // calculate the loop distance min
        double loop_distance_min( contact_type->GetPreferredDistanceRange().GetMin());

        // calculate the loop distance max
        // if the min is larger than max, which would happen if the loop is very short, then set it to min
        double loop_distance_max
        (
          std::max( std::min( contact_type->GetPreferredDistanceRange().GetMax(), max_loop_distance), loop_distance_min)
        );

        // calculate distance
        const double loop_distance( random::GetGlobalRandom().Random< double>( loop_distance_min, loop_distance_max));

        // length difference
        const double length_difference( ( NEIGHBOR_SSE.GetLength() - SELECTED_SSE.GetLength()) / 2.0);

        // calculate z_translation_amount
        const double z_translation_amount
        (
          // if selected SSE comes first the translation has to be -Z
          // if neighbor SSE is longer then translation should be positive
          double( math::ConvertBooleanToSign( !selected_sse_comes_first)) * length_difference
        );

        // calculate the twist angle
        const double tilt_angle( random::GetGlobalRandom().Double( contact_type->GetTiltAngleRange()));

        // find the translation to move the end of the SSE which will connect the loop to origin
        // if selected SSE comes first it's the C terminal edge, so we have to move down
        // otherwise it's the N terminal edge, so we have to move up
        linal::Vector3D translate_loop_end_to_center
        (
          0.0, 0.0, double( math::ConvertBooleanToSign( !selected_sse_comes_first)) * SELECTED_SSE.GetLength() / 2.0
        );

        // if not adding a strand
        if( SELECTED_SSE.GetType() != biol::GetSSTypes().STRAND)
        {
          // decide on a internal z rotation
          const double internal_z_rotation( random::GetGlobalRandom().Random< double>( 0.0, 2.0 * math::g_Pi));

          // rotate around Z
          transformation( math::RotationMatrix3D( coord::GetAxes().e_Z, internal_z_rotation));
        }

        // flip around flip axis
        transformation( math::RotationMatrix3D( flip_axis, math::g_Pi));

        // rotate the transformation by the tilt angle along the translation axis
        transformation( math::RotationMatrix3D( translation_direction, tilt_angle));

        // apply the translation
        transformation( translation_direction * loop_distance);

        // apply the z-translation
        transformation( linal::Vector3D( 0.0, 0.0, z_translation_amount));

        // apply the rotation of the neighbor SSE
        transformation( NEIGHBOR_SSE.GetOrientation());

        // return
        return storage::Pair< math::TransformationMatrix3D, bool>( transformation, true);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
