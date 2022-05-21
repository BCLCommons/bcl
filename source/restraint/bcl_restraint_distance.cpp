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
#include "restraint/bcl_restraint_distance.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Distance::s_Instance
    (
      GetObjectInstances().AddInstance( new Distance())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Distance::Distance() :
      m_Distance( util::GetUndefined< double>()),
      m_UpperBound( util::GetUndefined< double>()),
      m_LowerBound( util::GetUndefined< double>())
    {
    }

    //! @brief construct from a distance and a margin of error
    //! @param DISTANCE double which is the distance
    //! @param UPPER_BOUND is the upper bound possible for the distance
    //! @param LOWER_BOUND is the lower bound possible for the distance
    Distance::Distance( const double &DISTANCE, const double &UPPER_BOUND, const double &LOWER_BOUND) :
      m_Distance( DISTANCE),
      m_UpperBound( UPPER_BOUND),
      m_LowerBound( LOWER_BOUND)
    {
      BCL_Assert
      (
        m_UpperBound >= m_Distance && m_LowerBound <= m_Distance,
        "The lower bound " + util::Format()( LOWER_BOUND) + " should be lower than the distance " +
        util::Format()( DISTANCE) + " and the upper bound " + util::Format()( UPPER_BOUND) + " should be higher than the distance"
      );
    }

    //! @brief virtual copy constructor
    Distance *Distance::Clone() const
    {
      return new Distance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Distance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives formatted string describing this
    //! @return formatted string describing this
    const std::string Distance::GetIdentification() const
    {
      return util::Format()( LowerBound()) + " " +
        util::Format()( GetDistance()) + " " + util::Format()( UpperBound());
    }

    //! @brief GetDistance return m_Distance
    //! @return returns double which is m_Distance
    double Distance::GetDistance() const
    {
      return m_Distance;
    }

    //! @brief UpperBound returns m_UpperBound
    //! @return returns double which is m_UpperBound
    double Distance::UpperBound() const
    {
      return m_UpperBound;
    }

    //! @brief LowerBound returns m_LowerBound
    //! @return returns double which is m_LowerBound
    double Distance::LowerBound() const
    {
      return m_LowerBound;
    }

    //! @brief determines if the distance object is defined or not
    //! @return bool true if distance, upper bound, and lower bound are all defined
    bool Distance::IsDefined() const
    {
      return util::IsDefined( m_Distance) && util::IsDefined( m_LowerBound) && util::IsDefined( m_UpperBound);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read distance from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Distance::Read( std::istream &ISTREAM)
    {
      // write members
      io::Serialize::Read( m_Distance, ISTREAM);
      io::Serialize::Read( m_UpperBound, ISTREAM);
      io::Serialize::Read( m_LowerBound, ISTREAM);

      return ISTREAM;
    }

    //! @brief write distance to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &Distance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Distance, OSTREAM, INDENT);
      io::Serialize::Write( m_UpperBound, OSTREAM, INDENT);
      io::Serialize::Write( m_LowerBound, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates the upper error
    //! @return the upper error
    double Distance::GetUpperError() const
    {
      return m_UpperBound - m_Distance;
    }

    //! @brief calculates the lower error
    //! @return the lower error
    double Distance::GetLowerError() const
    {
      return m_Distance - m_LowerBound;
    }

  } // namespace restraint
} // namespace bcl

