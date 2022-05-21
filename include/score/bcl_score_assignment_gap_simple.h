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

#ifndef BCL_SCORE_ASSIGNMENT_GAP_SIMPLE_H_
#define BCL_SCORE_ASSIGNMENT_GAP_SIMPLE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AssignmentGapSimple
    //! @brief is for scoring a gap linearly as the gap size grows.
    //!
    //! @see @link example_score_assignment_gap_simple.cpp @endlink
    //! @author alexanns
    //! @date 04/23/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AssignmentGapSimple :
      public math::FunctionInterfaceSerializable< size_t, double>
    {

    private:

    //////////
    // data //
    //////////

      double m_Penalty; //!< penalty for gap

    public:

    //////////
    // data //
    //////////

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct ScoreGap from optional parameter
      AssignmentGapSimple( const double GAP_PENALTY);

      //! @brief Clone is the virtual Clone constructor
      //! @return a pointer to new AssignmentGapSimple which is a copy of this
      AssignmentGapSimple *Clone() const;

      //! destructor
      ~AssignmentGapSimple();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief access to penalty score value
      //! @return penalty score for one gap
      double GetPenalty() const;

      //! @brief change penalty score
      //! @param PENALTY the new penalty for a single gap
      void SetPenalty( const double &PENALTY);

      //! @brief evaluate the penalty for a gap with count SIZE
      //! @param SIZE size of the gap
      //! @return double which is the score for a gap of size SIZE
      double operator()( const size_t &SIZE) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class AssignmentGapSimple

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_ASSIGNMENT_GAP_SIMPLE_H_
