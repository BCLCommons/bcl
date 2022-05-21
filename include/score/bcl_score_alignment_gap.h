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

#ifndef BCL_SCORE_ALIGNMENT_GAP_H_
#define BCL_SCORE_ALIGNMENT_GAP_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_interface.h"
#include "align/bcl_align_sequence_interface.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignmentGap
    //! @brief is an alignment scoring function for scoring extension gaps
    //! This class scores the extension gaps of an alignment based on the length of alignment and its sequences.
    //! @tparam t_Member is the type of object that the Assignment will be storing
    //!
    //! @see @link example_score_alignment_gap.cpp @endlink
    //! @author heinzes1
    //! @date May 31, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignmentGap :
      public math::FunctionInterfaceSerializable< align::AlignmentInterface< t_Member>, double>
    {

    private:

    //////////
    // data //
    //////////

      double m_GapPenalty; //!< penalty for gap

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      //! @param GAP_PENALTY the penalty for a gap; default is 0.0
      AlignmentGap( const double GAP_PENALTY = 0.0) :
        m_GapPenalty( GAP_PENALTY)
      {
      }

      //! @brief Clone function
      //! @return pointer to new AlignmentGap
      AlignmentGap< t_Member> *Clone() const
      {
        return new AlignmentGap< t_Member>( *this);
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

      //! operator for scoring the alignment
      //! @param ALIGNMENT the alignment to be scored
      //! @return the gap extension score
      double operator()( const align::AlignmentInterface< t_Member> &ALIGNMENT) const
      {
        BCL_Exit( "AlignmentInterface<>::GetSequences is not implemented!", -1);
//        // initialize and count number of gaps for all sequences
//        size_t gaps( 0);
//        for
//        (
//          typename util::ShPtrList< align::SequenceInterface< t_Member> >::const_iterator
//            itr( ALIGNMENT.GetSequences().Begin()), itr_end( ALIGNMENT.GetSequences().End());
//          itr != itr_end;
//          ++itr
//        )
//        {
//          gaps += ( ALIGNMENT.GetSize() - ( *itr)->GetSize());
//        }
//
//        return m_GapPenalty * gaps; // calculate and return the gap score
        return 0;
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
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // template class AlignmentGap

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_ALIGNMENT_GAP_H_
