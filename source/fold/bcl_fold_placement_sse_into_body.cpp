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
#include "fold/bcl_fold_placement_sse_into_body.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "find/bcl_find_pick_criteria_interface.h"
#include "restraint/bcl_restraint_body.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PlacementSSEIntoBody::s_Instance
    (
      GetObjectInstances().AddInstance( new PlacementSSEIntoBody())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PlacementSSEIntoBody::PlacementSSEIntoBody() :
      m_BodyRestraint(),
      m_BodyPicker(),
      m_Orientation()
    {
    }

    //! @brief construct from the three member types
    //! @param RESTRAINT ShPtr to a restraint::Body which will be "m_BodyRestraint"
    //! @param BODY_PICKER ShPtr to a PickCriteriaInterface which will be "m_BodyPicker"
    //! @param ORIENTATION ShPtr to a FunctionInterface which will be "m_Orientation"
    PlacementSSEIntoBody::PlacementSSEIntoBody
    (
      const util::ShPtr< restraint::Body> &RESTRAINT,
      const util::ShPtr
      <
        find::PickCriteriaInterface< util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface>
      > &BODY_PICKER,
      util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > &ORIENTATION
    ) :
      m_BodyRestraint( RESTRAINT),
      m_BodyPicker( BODY_PICKER),
      m_Orientation( ORIENTATION)
    {
    }

    //! @brief virtual copy constructor
    PlacementSSEIntoBody *PlacementSSEIntoBody::Clone() const
    {
      return new PlacementSSEIntoBody( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PlacementSSEIntoBody::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetBodyRestraint returns a const reference to "m_BodyRestraint"
    //! @return returns a const reference to "m_BodyRestraint"
    const util::ShPtr< restraint::Body> &PlacementSSEIntoBody::GetBodyRestraint() const
    {
      return m_BodyRestraint;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Place returns transformatin matrix for the placement of an SSE into a ProteinModel
    //! @param SSE SSE which is to be placed
    //! @param PROTEIN_MODEL into which the SSE will be placed
    //! @return the transformationmatrix3d to place SSE in PROTEIN_MODEL
    storage::Pair< math::TransformationMatrix3D, bool> PlacementSSEIntoBody::Place
    (
      const assemble::SSE &SSE,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // make sure that "m_BodyPicker" is defined
      BCL_Assert( m_BodyPicker.IsDefined(), "m_BodyPicker ShPtr is not defined");

      // return storage::Pair given by DetermineTransformatinMatrix3DandBool given a body picked by "m_BodyPicker"
      // out of the bodies that are not occupied by an of the SSEs in "PROTEIN_MODEL"
      return DetermineTransformatinMatrix3DandBool
      (
        m_BodyPicker->Pick( m_BodyRestraint->GetUnoccupied( PROTEIN_MODEL.GetSSEs()), SSE)
      );
    }

    //! returns placement for the t_ObjectType without taking ProteinModel into consideration i.e. if it is empty
    //! @param SSE SSE which is to be placed
    //! @return the transformationmatrix3d to place SSE in PROTEIN_MODEL
    storage::Pair< math::TransformationMatrix3D, bool> PlacementSSEIntoBody::Place
    (
      const assemble::SSE &SSE
    ) const
    {
      // make sure that "m_BodyPicker" is defined
      BCL_Assert( m_BodyPicker.IsDefined(), "m_BodyPicker ShPtr is not defined");

      // return storage::Pair given by DetermineTransformatinMatrix3DandBool given a body picked by "m_BodyPicker"
      // and "SSE"
      return DetermineTransformatinMatrix3DandBool( m_BodyPicker->Pick( *m_BodyRestraint->GetBody(), SSE));
    }

    //! @brief DetermineTransformationMatrix3DandBool figures out what the transformation matrix and bool should be
    //! @param RESTRAINT_BODY that was chosen to be the body into which the SSE will be placed
    //! @return the transformationmatrix3d to place SSE in PROTEIN_MODEL and a bool whether successful or not
    storage::Pair< math::TransformationMatrix3D, bool> PlacementSSEIntoBody::DetermineTransformatinMatrix3DandBool
    (
      const util::ShPtr< assemble::SSEGeometryInterface> &RESTRAINT_BODY
    ) const
    {
      // true if "RESTRAINT_BODY" is not defined meaning no "assemble::SSEGeometryInterface" could be picked by "m_BodyRestraint"
      if( !RESTRAINT_BODY.IsDefined())
      {
        // write message
        BCL_MessageDbg( "RESTRAINT_BODY is not defined\n");

        // no "assemble::SSEGeometryInterface" could be picked by "m_BodyRestraint" so return empty math::TransformationMatrix3D()
        // and a boolean of false
        return storage::Pair< math::TransformationMatrix3D, bool>( math::TransformationMatrix3D(), false);
      }
      else //< "RESTRAINT_BODY" is defined so there is a restraint body to be placed into
      {
        // write message
        BCL_MessageDbg
        (
          "RESTRAINT_BODY is defined\n"
        );

        // output body extents of the "RESTRAINT_BODY"
        BCL_MessageDbg
        (
          "body is defined with main axis: \n " +
          util::Format()( RESTRAINT_BODY->GetMainAxis())
        );

        // create storage::Pair "result_pair" and initialize with TransformationMatrix3D for putting an SSE into the
        // restraint body and true bool
        storage::Pair< math::TransformationMatrix3D, bool> result_pair
        (
          // the TransformationMatrix3D is determined by "m_Orientation" according to "RESTRAINT_BODY"
          math::TransformationMatrix3D( *m_Orientation->operator()( RESTRAINT_BODY->GetOrientation()).GetArgument()),
          true
        );

        // end
        return result_pair;
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PlacementSSEIntoBody::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_BodyRestraint, ISTREAM);
      io::Serialize::Read( m_BodyPicker, ISTREAM);
      io::Serialize::Read( m_Orientation, ISTREAM);

      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @return ostream which was read from
    std::ostream &PlacementSSEIntoBody::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // read members
      io::Serialize::Write( m_BodyRestraint, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BodyPicker, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Orientation, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
