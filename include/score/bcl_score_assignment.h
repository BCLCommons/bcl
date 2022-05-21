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

#ifndef BCL_SCORE_ASSIGNMENT_H_
#define BCL_SCORE_ASSIGNMENT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "align/bcl_align_assignment.h"
#include "function/bcl_function_binary_interface.h"
#include "function/bcl_function_unary_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Assignment
    //! @brief scores an assignment of t_Members
    //! This class scores an assignment by using MemberPairScores (FunctionInterface with two t_Members and double)
    //! @tparam t_Member is the type of object that the Assignment stores
    //!
    //! @see @link example_score_assignment.cpp @endlink
    //! @author heinzes1, woetzen
    //! @date May 28, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class Assignment :
      public function::UnaryInterface< const align::Assignment< t_Member>, double>
    {

    private:

    //////////
    // data //
    //////////

      // internal typedef
      typedef function::BinaryInterface< const t_Member, const t_Member, double> MemberPairScore;

      //! Function to score pair of members
      util::ShPtr< MemberPairScore> m_MemberPairScore;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      Assignment() :
        m_MemberPairScore()
      {
      }

      //! @brief constructor taking a MemberPairScore
      //! @param SCORE_PAIR ShPtr to member pair score
      Assignment( const util::ShPtr< MemberPairScore> &SCORE_PAIR) :
        m_MemberPairScore( SCORE_PAIR)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Assignment
      Assignment< t_Member> *Clone() const
      {
        return new Assignment< t_Member>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ASSIGNMENT and returning its score
      //! @param ASSIGNMENT to be evaluated
      //! @return assignment score
      double operator()( const align::Assignment< t_Member> &ASSIGNMENT) const
      {
        // check for at least two members; normalization does not work otherwise
        const util::SiPtrList< const t_Member> members( ASSIGNMENT.GetMembers()); // initialize
        if( members.GetSize() < 2)
        {
          return 0.0;
        }

        double score( 0.0); // initialize
        for // loop over all pairs of members
        (
          typename util::SiPtrList< const t_Member>::const_iterator itr_a( members.Begin()), itr_end( members.End());
          itr_a != itr_end;
          ++itr_a
        )
        {
          if( !itr_a->IsDefined())
          {
            continue;
          }

          const t_Member &member_a( **itr_a);

          typename util::SiPtrList< const t_Member>::const_iterator itr_b( itr_a); // make itr_b a copy of itr_a
          for( ++itr_b; itr_b != itr_end; ++itr_b) // increase itr_b by one to start after itr_a
          {
            if( !itr_b->IsDefined())
            {
              continue;
            }
            const t_Member &member_b( **itr_b);

            score += m_MemberPairScore->operator()( member_a, member_b);
          }
        }

        // return normalized score
        return score / double( members.GetSize() * ( members.GetSize() - 1) / 2);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_MemberPairScore, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_MemberPairScore, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class Assignment

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_ASSIGNMENT_H_ 
