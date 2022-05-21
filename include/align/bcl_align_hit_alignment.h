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

#ifndef BCL_ALIGN_HIT_ALIGNMENT_H_
#define BCL_ALIGN_HIT_ALIGNMENT_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HitAlignment
    //! @brief TODO: add a brief comment
    //! @details TODO: add an general comment to this class
    //!
    //! @see @link example_align_hit_alignment.cpp @endlink
    //! @author heinzes1
    //! @date Jan 12, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class HitAlignment :
//      public AlignmentSimple< t_Member>
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

//        // the two parent alignments that the hit combines
//        util::SiPtr< const AlignmentSimple< t_Member> > m_ParentAlignmentA;
//        util::SiPtr< const AlignmentSimple< t_Member> > m_ParentAlignmentB;

      // begin and end position of hit alignment in the two parent aligments; used for finding a nonoverlapping set
      size_t m_BeginA;
      size_t m_EndA;
      size_t m_BeginB;
      size_t m_EndB;

      // the score of the hit alignment using the score::Assignment; used for applying a threshold
      double m_Score;

    public:

    //////////
    // data //
    //////////

//      //! single instance of that class
//      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //
      HitAlignment() :
//        Alignment< t_Member>(),
//        m_ParentAlignmentA(),
//        m_ParentAlignmentB(),
        m_BeginA( 0),
        m_EndA( 0),
        m_BeginB( 0),
        m_EndB( 0),
        m_Score( 0.0)
      {
      }

//      //
//      template< typename t_Iterator>
//      HitAlignment
//      (
//        const t_Iterator &FIRST,
//        const t_Iterator &LAST,
//        const util::SiPtr< const AlignmentInterface< t_Member> > &PARENT_ALIGNMENT,
//        const size_t BEGIN,
//        const size_t END
//      ) :
////        Alignment< t_Member>( FIRST, LAST),
////        m_ParentAlignmentA( PARENT_ALIGNMENT),
////        m_ParentAlignmentB( PARENT_ALIGNMENT),
//        m_BeginA( BEGIN),
//        m_EndA( END),
//        m_BeginB( BEGIN),
//        m_EndB( END),
//        m_Score( 0.0)
//      {
//      }

//      //
//      template< typename t_Iterator>
//      HitAlignment
//      (
//        const t_Iterator &FIRST_A,
//        const t_Iterator &LAST_A,
//        const util::SiPtr< const AlignmentInterface< t_Member> > &PARENT_ALIGNMENT_A,
//        const size_t BEGIN_A,
//        const size_t END_A,
//        const t_Iterator &FIRST_B,
//        const t_Iterator &LAST_B,
//        const util::SiPtr< const AlignmentInterface< t_Member> > &PARENT_ALIGNMENT_B,
//        const size_t BEGIN_B,
//        const size_t END_B
//      ) :
////        Alignment< t_Member>(),
////        m_ParentAlignmentA( PARENT_ALIGNMENT_A),
////        m_ParentAlignmentB( PARENT_ALIGNMENT_B),
//        m_BeginA( BEGIN_A),
//        m_EndA( END_A),
//        m_BeginB( BEGIN_B),
//        m_EndB( END_B),
//        m_Score( 0.0)
//      {
//        // create new shptr to assignments, insert t_Members and insert shptr into alignment (base class)
//        for( t_Iterator itr_a( FIRST_A), itr_b( FIRST_B); itr_a != LAST_A; ++itr_a, ++itr_b)
//        {
//          util::ShPtr< Assignment< t_Member> > ass( new Assignment< t_Member>());
//          ass->GetMembers().InsertElements( ass->GetMembers().End(), ( *itr_a)->GetMembers());
//          ass->GetMembers().InsertElements( ass->GetMembers().End(), ( *itr_b)->GetMembers());
//          AlignmentSimple< t_Member>::PushBack( ass);
//        }
//
//        //debug
//        BCL_MessageStd( "DEBUG HitAlignment::HitAlignment()");
//        for
//        (
//          typename util::ShPtrList< Assignment< t_Member> >::const_iterator itr( util::ShPtrList< Assignment< t_Member> >::Begin());
//          itr != util::ShPtrList< Assignment< t_Member> >::End();
//          ++itr
//        )
//        {
//          util::ShPtr< Assignment< t_Member> > x2( *itr);
//          for
//          (
//            typename util::SiPtrList< const t_Member>::iterator itr_member( ( *x2).GetMembers().Begin());
//            itr_member != ( *x2).GetMembers().End();
//            ++itr_member
//          )
//          {
//            util::SiPtr< const t_Member> member_ptr( *itr_member);
//            util::SiPtr< const biol::AABase> aabase_ptr( member_ptr);
//            std::string str( aabase_ptr->GetIdentification());
//            BCL_MessageStd( "DEBUG aabase=" + util::Format()( str));
//          }
//        }
//        BCL_MessageStd( "DEBUG m_ParentAlignmentA=" + util::Format()( &m_ParentAlignmentA));
//        BCL_MessageStd( "DEBUG m_ParentAlignmentB=" + util::Format()( &m_ParentAlignmentB));
//        BCL_MessageStd( "DEBUG m_Score=" + util::Format()( m_Score) + "|m_BeginA=" + util::Format()( m_BeginA) + "|m_EndA=" + util::Format()( m_EndA) + "|m_BeginB=" + util::Format()( m_BeginB) + "|m_EndB=" + util::Format()( m_EndB));
//      }

      //! @brief Clone function
      //! @return pointer to new Hit
      HitAlignment *Clone() const
      {
        return new HitAlignment< t_Member>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

//      //
//      util::SiPtr< const AlignmentSimple< t_Member> > GetParentAlignmentA() const
//      {
//        return m_ParentAlignmentA;
//      }
//
//      //
//      util::SiPtr< const AlignmentSimple< t_Member> > GetParentAlignmentB() const
//      {
//        return m_ParentAlignmentB;
//      }

      //
      size_t GetBeginInParentAlignmentA() const
      {
        return m_BeginA;
      }

      //
      size_t GetEndInParentAlignmentA() const
      {
        return m_EndA;
      }

      //
      size_t GetBeginInParentAlignmentB() const
      {
        return m_BeginB;
      }

      //
      size_t GetEndInParentAlignmentB() const
      {
        return m_EndB;
      }

      //
      double GetScore() const
      {
        return m_Score;
      }

      //
      void SetScore( double SCORE)
      {
        m_Score = SCORE;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //
      bool operator<( const HitAlignment &HIT) const
      {
        return m_Score < HIT.m_Score;
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
        // exit if function is called
        BCL_Exit( "can't read in Assignments to create Alignment", -1);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
//        // write base class
//        AlignmentSimple< t_Member>::Write( OSTREAM, INDENT);

        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class HitAlignment
  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_HIT_ALIGNMENT_H_
