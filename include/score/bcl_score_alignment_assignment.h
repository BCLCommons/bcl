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

#ifndef BCL_SCORE_ALIGNMENT_ASSIGNMENT_H_
#define BCL_SCORE_ALIGNMENT_ASSIGNMENT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_interface.h"
#include "function/bcl_function_unary_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignmentAssignment
    //! @brief scores all assignment in an alignment (vertical score)
    //! @detail scores all assignments in an alignment and sums up all the score
    //! @tparam t_Member type of elements aligned
    //!
    //! @see @link example_score_alignment_assignment.cpp @endlink
    //! @author heinzes1, woetzen
    //! @date May 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignmentAssignment :
      public function::UnaryInterface< const align::AlignmentInterface< t_Member>, double>
    {

    private:

    //////////
    // data //
    //////////

      //! the assignment scoring function
      util::ShPtr< function::UnaryInterface< const align::Assignment< t_Member>, double> > m_AssignmentScore;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking an assignment score
      //! @param SP_ASSIGNMENT_SCORE ShPtr to an assignment score
      AlignmentAssignment
      (
        const util::ShPtr< function::UnaryInterface< const align::Assignment< t_Member>, double> > &SP_ASSIGNMENT_SCORE
      ) :
        m_AssignmentScore( SP_ASSIGNMENT_SCORE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new AlignmentVertically< t_Member>
      AlignmentAssignment< t_Member> *Clone() const
      {
        return new AlignmentAssignment< t_Member>( *this);
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

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! score alignment
      //! @param ALIGNMENT the alignment to be scored
      //! @return the score
      double operator()( const align::AlignmentInterface< t_Member> &ALIGNMENT) const
      {
        double score( 0);

        for
        (
          typename align::AlignmentInterface< t_Member>::const_iterator itr( ALIGNMENT.GetAssignments().Begin()), itr_end( ALIGNMENT.GetAssignments().End());
          itr != itr_end;
          ++itr
        )
        {
          score += m_AssignmentScore->operator()( **itr);
        }

        return score;
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
        // read members
        io::Serialize::Read( m_AssignmentScore, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_AssignmentScore, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class AlignmentAssignment

    // instantiate s_Instance
    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> AlignmentAssignment< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignmentAssignment< t_Member>())
    );

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_ALIGNMENT_ASSIGNMENT_H_ 
