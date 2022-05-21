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
#include "fold/bcl_fold_mutate_protein_model_sse_swap_body.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSESwapBody::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelSSESwapBody())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSESwapBody::MutateProteinModelSSESwapBody() :
      m_BodyRestraint(),
      m_BodyPicker(),
      m_Locator(),
      m_Scheme( GetStaticClassName< MutateProteinModelSSESwapBody>())
    {
    }

    //! @brief constructor taking each of the member variable types
    //! @param BODY_RESTRAINT ShPtr to a restraint::Body which will be "m_BodyRestraint"
    //! @param BODY_PICKER ShPtr to a PickerInterface which will be "m_BodyPicker"
    //! @param LOCATOR ShPtr to a LocatorInterface which will by "m_Locator"
    MutateProteinModelSSESwapBody::MutateProteinModelSSESwapBody
    (
      const util::ShPtr< restraint::Body> &BODY_RESTRAINT,
      const util::ShPtr
      <
        find::PickCriteriaInterface< util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface>
      > &BODY_PICKER,
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &LOCATOR,
      const std::string &SCHEME
    ) :
      m_BodyRestraint( BODY_RESTRAINT),
      m_BodyPicker( BODY_PICKER),
      m_Locator( LOCATOR),
      m_Scheme( SCHEME)
    {
        BCL_Assert( m_BodyPicker.IsDefined(), "body picker pointer that is passed in constructor needs to be defined for move to work");
    }

    //! @brief clone
    MutateProteinModelSSESwapBody *MutateProteinModelSSESwapBody::Clone() const
    {
      return new MutateProteinModelSSESwapBody( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSESwapBody::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel>
    MutateProteinModelSSESwapBody::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // true if PROTEIN_MODEL has no sses in it
      if( PROTEIN_MODEL.GetSSEs().IsEmpty())
      {
        // then return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // create ShPtr to an SSE "new_sse" and initialize with a Clone of the SSE located by "m_Locator" out of
      // PROTEIN_MODEL. "new_sse" will be swapped from its current position into a restraint body
      util::ShPtr< assemble::SSE> new_sse( m_Locator->Locate( PROTEIN_MODEL)->Clone());

      BCL_Assert( m_BodyPicker.IsDefined(), "body picker pointer needs to be defined for move to work");

      // create ShPtr to a assemble::SSEGeometryInterface "unoccupied_body" and initialize with the assemble::SSEGeometryInterface chosen by "m_BodyPicker"
      // out of the bodies that are not occupied by any of the SSE's in PROTEIN_MODEL
      // this is the body which "new_sse" will be moved into
      const util::ShPtr< assemble::SSEGeometryInterface> unoccupied_body
      (
        m_BodyPicker->Pick( m_BodyRestraint->GetUnoccupied( PROTEIN_MODEL.GetSSEs()), *new_sse)
      );

      // true if unoccupied_body is not defined (i.e. all density rods are occupied by sses of the protein model,
      // there is no empty density rod)
      if( !unoccupied_body.IsDefined())
      {
        // then return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // create TransformationMatrix3D "transform_new_sse" and initialize with the inverse transformation matrix
      // of "new_sse" which "moves" the transformation matrix to the origin so
      // rotations of "new_sse" will occur internally
      math::TransformationMatrix3D transform_new_sse( math::Inverse( new_sse->GetOrientation()));

      // apply the transformation matrix of "unoccupied_body" to "transform_new_sse"
      transform_new_sse( unoccupied_body->GetOrientation());

      // move "new_sse" to the position of "unoccupied_body" based on "transform_new_sse"
      new_sse->Transform( transform_new_sse);

      // create ShPtr to a ProteinModel "new_model" and initialize with a clone of "PROTEIN_MODEL"
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // replace "sse" in "new_model" with "new_sse"
      new_model->Replace( new_sse);

      // return MutateResult constructed from "new_model" and this mutate object
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelSSESwapBody::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_BodyRestraint, ISTREAM);
      io::Serialize::Read( m_BodyPicker, ISTREAM);
      io::Serialize::Read( m_Locator, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &MutateProteinModelSSESwapBody::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_BodyRestraint, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BodyPicker, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Locator, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
