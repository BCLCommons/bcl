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
#include "fold/bcl_fold_placement_sse_next_to_sse.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "contact/bcl_contact_types.h"
#include "find/bcl_find_locator_criteria_interface.h"
#include "io/bcl_io_serialization.h"
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
    const util::SiPtr< const util::ObjectInterface> PlacementSSENextToSSE::s_Instance
    (
      util::Enumerated< PlacementInterface< assemble::SSE, assemble::ProteinModel> >::AddInstance( new PlacementSSENextToSSE())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PlacementSSENextToSSE::PlacementSSENextToSSE() :
      m_NeighborSSELocator()
    {
    }

    //! @brief constructor from a LocatorCriteriaInterface
    //! @param LOCATOR_CRITERIA LocatorCriteriaInterface to be used
    PlacementSSENextToSSE::PlacementSSENextToSSE
    (
      const find::LocatorCriteriaInterface
      <
        util::SiPtr< const assemble::SSE>,
        assemble::DomainInterface,
        assemble::SSE
      > &LOCATOR_CRITERIA
    ) :
      m_NeighborSSELocator( LOCATOR_CRITERIA.Clone())
    {
    }

    //! @brief constructor from a ShPtr to LocatorCriteriaInterface
    //! @param SP_LOCATOR_CRITERIA ShPtr to LocatorCriteriaInterface to be used
    PlacementSSENextToSSE::PlacementSSENextToSSE
    (
      const util::ShPtr
      <
        find::LocatorCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>,
          assemble::DomainInterface,
          assemble::SSE
        >
      > &SP_LOCATOR_CRITERIA
    ) :
      m_NeighborSSELocator( *SP_LOCATOR_CRITERIA)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PlacementSSENextToSSE
    PlacementSSENextToSSE *PlacementSSENextToSSE::Clone() const
    {
      return new PlacementSSENextToSSE( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PlacementSSENextToSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &PlacementSSENextToSSE::GetAlias() const
    {
      static const std::string s_alias( "PlacementSSENextToSSE");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PlacementSSENextToSSE::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Places an SSE next to other SSEs in a protein model.");
      serializer.AddInitializer
      (
        "locator",
        "locator for neighboring SSEs",
        io::Serialization::GetAgent( &m_NeighborSSELocator)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate placement for the given SSE at a random orientation wrt to a located SSE from the given model
    //! @param SELECTED_SSE SiPtr to SSE to be placed, it needs to be in the origin
    //! @param PROTEIN_MODEL to which the SSE is going to be added
    storage::Pair< math::TransformationMatrix3D, bool> PlacementSSENextToSSE::Place
    (
      const assemble::SSE &SELECTED_SSE,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize neighbor_sse
      util::SiPtr< const assemble::SSE> neighbor_sse;

      // if there is only one sse in the model set it to neighbor_sse
      if( PROTEIN_MODEL.GetSSEs().GetSize() == 1)
      {
        neighbor_sse = PROTEIN_MODEL.GetSSEs().FirstElement();
      }
      // else if more than one SSE
      else
      {
        // locate the neighbor sse using m_NeighborSSELocator
        neighbor_sse = m_NeighborSSELocator->Locate( PROTEIN_MODEL, SELECTED_SSE);
      }

      // if the neighbor sse is empty
      if( !neighbor_sse.IsDefined())
      {
        // give a warning and return
        BCL_MessageCrt( "No neighbor sses are located!");
        return storage::Pair< math::TransformationMatrix3D, bool>( math::TransformationMatrix3D(), false);
      }

      BCL_MessageVrb( "neighbor SSE selected: " + neighbor_sse->GetIdentification());

      // calculate the transformation and return
      return Place( SELECTED_SSE, *neighbor_sse);
    }

    //! @brief generate placement for the given SSE at a random orientation wrt to a located SSE from the given model
    //! @param SELECTED_SSE SiPtr to SSE to be placed
    //! @param NEIGHBOR_SSE reference SSE
    storage::Pair< math::TransformationMatrix3D, bool> PlacementSSENextToSSE::Place
    (
      const assemble::SSE &SELECTED_SSE,
      const assemble::SSE &NEIGHBOR_SSE
    )
    {
      // initialize the transformation matrix that will be applied to this sse
      math::TransformationMatrix3D transformation;

      // decide whether to flip
      if( random::GetGlobalRandom().Boolean())
      {
        // apply a flip around a random axis
        transformation( coord::Axis( random::GetGlobalRandom().Random< size_t>( 2)), math::g_Pi);
      }

      // decide positive or negative translation along random_axis
      const int translation_sign( random::GetGlobalRandom().Sign());

      // create the contact type
      contact::Type contact_type( contact::GetTypes().TypeFromSSTypes( SELECTED_SSE, NEIGHBOR_SSE));

      // initialize translation vector and the distance
      linal::Vector3D translation_direction;

      // if helix helix
      if( contact_type == contact::GetTypes().HELIX_HELIX)
      {
        // determine a random angle between -90 and +90
        const double angle( random::GetGlobalRandom().Random< double>( -math::g_Pi / 2.0, math::g_Pi / 2.0));

        // use this angle to determine a random unit vector on xy plane and scale the translation
        translation_direction = linal::Vector3D( cos( angle), sin( angle), 0.0);
      }
      // if the selected_sse is a helix and the neighbor sse is strand
      else if( contact_type == contact::GetTypes().HELIX_SHEET)
      {
        // then we can only use x axis
        translation_direction = linal::Vector3D( 1.0, 0.0, 0.0);
      }
      // if the selected SSE is a strand and the neighbor sse is helix
      else if( contact_type == contact::GetTypes().SHEET_HELIX)
      {
        // determine a random angle between -90 and +90
        const double angle( random::GetGlobalRandom().Random< double>( -math::g_Pi / 2.0, math::g_Pi / 2.0));

        // use this angle to determine a random unit vector on xy plane
        translation_direction = linal::Vector3D( cos( angle), sin( angle), 0.0);

        // we also need a rotation around Z to make sure the strand side chains are facing the helix
        transformation( math::RotationMatrix3D( coord::GetAxes().e_Z, -angle));

      }
      else if( contact_type == contact::GetTypes().STRAND_STRAND)
      {
        // determine randomly to use x (SHEET_SHEET) or y axis( STRAND_STRAND)
        const coord::Axis
            trans_axis( random::GetGlobalRandom().Boolean() ? coord::GetAxes().e_X : coord::GetAxes().e_Y);

        // set the translation direction
        translation_direction( trans_axis) = 1.0;

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
        BCL_MessageCrt
        (
          "undefined contact type! " + contact_type.GetName() +
          " selected SSE: " + SELECTED_SSE.GetIdentification() +
          " neighbor SSE: " + NEIGHBOR_SSE.GetIdentification()
        );
        return storage::Pair< math::TransformationMatrix3D, bool>( math::TransformationMatrix3D(), false);
      }

      // get a random distance within the range specified by the contact type
      const double distance( random::GetGlobalRandom().Double( contact_type->GetPreferredDistanceRange()));

      // get a random tilt angle within the range allowed by the contact type
      const double tilt_angle( random::GetGlobalRandom().Double( contact_type->GetTiltAngleRange()));

      // decide positive or negative internal translation along z axis
      const int z_translation_sign( random::GetGlobalRandom().Sign());

      // length difference
      const double length_difference( NEIGHBOR_SSE.GetLength() - SELECTED_SSE.GetLength());

      // initialize internal z translation
      linal::Vector3D z_translation
      (
        0.0, 0.0,
        double( z_translation_sign) * random::GetGlobalRandom().Double
        (
          math::Range< double>
          (
            std::min( 0.0, length_difference / 2.0),
            std::max( length_difference / 2.0, 0.0)// half the length difference
          )
        )
      );

      // if not adding a strand
      if( SELECTED_SSE.GetType() != biol::GetSSTypes().STRAND)
      {
        // decide on a internal z rotation
        const double internal_z_rotation( random::GetGlobalRandom().Random< double>( 0.0, 2 * math::g_Pi));

        // rotate around Z
        transformation( math::RotationMatrix3D( coord::GetAxes().e_Z, internal_z_rotation));
      }

      // rotate the transformation by the tilt angle along the translation axis
      transformation( math::RotationMatrix3D( translation_direction, tilt_angle));

      // form the translation vector
      linal::Vector3D translation( translation_direction * distance * translation_sign);

      // now apply the internal z translation
      transformation( z_translation);

      // apply the calculated translation to transformation
      transformation( translation);

      // apply the orientation of neighbor sse to transformation
      transformation( NEIGHBOR_SSE.GetOrientation());

      // return
      return storage::Pair< math::TransformationMatrix3D, bool>( transformation, true);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
