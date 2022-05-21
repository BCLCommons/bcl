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

#ifndef BCL_RESTRAINT_DISTANCE_H_
#define BCL_RESTRAINT_DISTANCE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Distance
    //! @brief This is a class for holding distance data. It represents a distance restraint with upper and lower bounds
    //!
    //! @see @link example_restraint_distance.cpp @endlink
    //! @author alexanns
    //! @date 05/12/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API Distance :
      public util::ObjectInterface
    {
    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      double m_Distance;   //!< the distance
      double m_UpperBound; //!< the upper bound possible for the distance
      double m_LowerBound; //!< the lower bound possible for the distance

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Distance();

      //! @brief construct from a distance and a margin of error
      //! @param DISTANCE double which is the distance
      //! @param UPPER_BOUND is the upper bound possible for the distance
      //! @param LOWER_BOUND is the lower bound possible for the distance
      Distance( const double &DISTANCE, const double &UPPER_BOUND, const double &LOWER_BOUND);

      //! @brief virtual copy constructor
      Distance *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives formatted string describing this
      //! @return formatted string describing this
      const std::string GetIdentification() const;

      //! @brief GetDistance return m_Distance
      //! @return returns double which is m_Distance
      double GetDistance() const;

      //! @brief UpperBound returns m_UpperBound
      //! @return returns double which is m_UpperBound
      double UpperBound() const;

      //! @brief LowerBound returns m_LowerBound
      //! @return returns double which is m_LowerBound
      double LowerBound() const;

      //! @brief determines if the distance object is defined or not
      //! @return bool true if distance, upper bound, and lower bound are all defined
      bool IsDefined() const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read distance from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write distance to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////
    public:

      //! @brief calculates the upper error
      //! @return the upper error
      double GetUpperError() const;

      //! @brief calculates the lower error
      //! @return the lower error
      double GetLowerError() const;
    };

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_DISTANCE_H_
