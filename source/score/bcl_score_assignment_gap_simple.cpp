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
#include "score/bcl_score_assignment_gap_simple.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AssignmentGapSimple::s_Instance
    (
      GetObjectInstances().AddInstance( new AssignmentGapSimple( 0.0))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor from optional parameter
    AssignmentGapSimple::AssignmentGapSimple( const double GAP_PENALTY) :
      m_Penalty( GAP_PENALTY)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new AssignmentGapSimple which is a copy of this
    AssignmentGapSimple *AssignmentGapSimple::Clone() const
    {
      return new AssignmentGapSimple( *this);
    }

    //! destructor
    AssignmentGapSimple::~AssignmentGapSimple()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AssignmentGapSimple::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function GetData returns "m_Penalty"
    //! @return "m_Penalty"
    double AssignmentGapSimple::GetPenalty() const
    {
      return m_Penalty;
    }

    //! @brief function SetData changes "m_Penalty" to a new penalty
    void AssignmentGapSimple::SetPenalty( const double &PENALTY)
    {
      m_Penalty = PENALTY;
    }

    //! @brief operator () overwritten for the FunctionInterface
    //! @param SIZE size_t which is the size of the gap
    //! @return double which is the score for a gap of size SIZE
    double AssignmentGapSimple::operator()( const size_t &SIZE) const
    {
      // return "m_Penalty" multiplied by SIZE
      return m_Penalty * SIZE;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AssignmentGapSimple::Read( std::istream &ISTREAM)
    {
      // read in member
      io::Serialize::Read( m_Penalty, ISTREAM);

      // return ISTREAM
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &AssignmentGapSimple::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Penalty, OSTREAM, INDENT) << '\n';

      // return OSTREAM
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl

