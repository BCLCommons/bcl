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
#include "fold/bcl_fold_mutate_sheet_twist.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "io/bcl_io_serialization.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateSheetTwist::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSheetTwist())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateSheetTwist::MutateSheetTwist() :
      m_AngleRange( math::Range< double>( 0.0, math::Angle::Radian( 10.0)))
    {
    }

    //! @brief constructor from an angle range
    //! @param ANGLE_RANGE range of changes in the twist angles in radians to be applied
    MutateSheetTwist::MutateSheetTwist( const math::Range< double> ANGLE_RANGE) :
      m_AngleRange( ANGLE_RANGE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateSheetTwist
    MutateSheetTwist *MutateSheetTwist::Clone() const
    {
      return new MutateSheetTwist( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateSheetTwist::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateSheetTwist::GetAlias() const
    {
      static const std::string s_name( "MutateSheetTwist");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateSheetTwist::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Adjusts twist angles in a sheet to random values.");
      serializer.AddInitializer
      (
        "angle range",
        "range of allowed twist angles",
        io::Serialization::GetAgent( &m_AngleRange)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes a Sheet and return a mutated Sheet
    //! @param SHEET Sheet which will be mutated
    //! @return MutateResult with the mutated Sheet
    math::MutateResult< assemble::Domain> MutateSheetTwist::operator()( const assemble::Domain &SHEET) const
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

      // iterate over SSEs in this sheet using the SSE vector
      for
      (
        util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
          sse_itr( SHEET.GetTopology()->GetElements().Begin()), sse_itr_end( SHEET.GetTopology()->GetElements().End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // get the next strand
        const assemble::SSEGeometryInterface &this_strand( **sse_itr);

        // pick randomly the angle change to be applied
        const double angle_change
        (
          random::GetGlobalRandom().Double( m_AngleRange) * double( random::GetGlobalRandom().Sign())
        );
        BCL_MessageVrb( "Angle selected " + util::Format()( angle_change));

        // construct the transformation need for the angle change
        // first step is to move to the origin
        math::TransformationMatrix3D transformation( math::Inverse( this_strand.GetOrientation()));
        // rotate with the angle change along the y axis
        transformation( coord::GetAxes().e_Y, angle_change);
        // move back to the original location
        transformation( this_strand.GetOrientation());

        // make a copy of this strand and apply the transformation
        // since this is a dynamic we have to make sure it was valid
        util::ShPtr< assemble::SSE> new_strand( this_strand.Clone());
        BCL_Assert( new_strand.IsDefined(), "The dynamic cast from SSEGeometryInterface to SSE failed!");
        new_strand->Transform( transformation);

        // replace the corresponding strand with the new strand in the new sheet
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
    std::istream &MutateSheetTwist::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AngleRange, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateSheetTwist::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // read members
      io::Serialize::Write( m_AngleRange, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
