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

#ifndef BCL_ALIGN_ASSIGNMENT_H_
#define BCL_ALIGN_ASSIGNMENT_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Assignment
    //! @brief class stores objects assigned together
    //! @details This is a template class for storing objects which go together/are assigned together for some reason.
    //!
    //! @see @link example_align_assignment.cpp @endlink
    //! @author heinzes1
    //! @date 10/12/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class Assignment :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      util::SiPtrList< const t_Member> m_Members; //!< m_Members holds the t_Members

      //! default gap char for empty SiPtr
      static const char s_DefaultGapChar = '-';

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Assignment()
      {
      }

      //! @brief construct assignment from a member
      //! @param SP_MEMBER first member of assignment
      Assignment( const util::SiPtr< const t_Member> &SP_MEMBER) :
        m_Members( 1, SP_MEMBER)
      {
      }

      //! @brief construct assignment from two members
      //! @param SP_MEMBER_A first member of assignment
      //! @param SP_MEMBER_B second member of assignment
      Assignment( const util::SiPtr< const t_Member> &SP_MEMBER_A, const util::SiPtr< const t_Member> &SP_MEMBER_B) :
        m_Members()
      {
        m_Members.PushBack( SP_MEMBER_A);
        m_Members.PushBack( SP_MEMBER_B);
      }

      //! @brief construct from two assignments
      //! @param ASSIGNMENT_TOP top assignment
      //! @param ASSIGNMENT_BOTTOM bottom assignment
      Assignment( const Assignment< t_Member> &ASSIGNMENT_TOP, const Assignment< t_Member> &ASSIGNMENT_BOTTOM) :
        m_Members( ASSIGNMENT_TOP.m_Members)
      {
        m_Members.Append( ASSIGNMENT_BOTTOM.m_Members);
      }

      //! @brief construct from member and assignment
      //! @param SP_MEMBER first member of assignment
      //! @param ASSIGNMENT assignments at bottom
      Assignment( const util::SiPtr< const t_Member> &SP_MEMBER, const Assignment< t_Member> &ASSIGNMENT) :
        m_Members( ASSIGNMENT.m_Members)
      {
        m_Members.PushFront( SP_MEMBER);
      }

      //! @brief construct from assignment and member
      //! @param ASSIGNMENT assignments at top
      //! @param SP_MEMBER bottom member of assignment
      Assignment( const Assignment< t_Member> &ASSIGNMENT, const util::SiPtr< const t_Member> &SP_MEMBER) :
        m_Members( ASSIGNMENT.m_Members)
      {
        m_Members.PushBack( SP_MEMBER);
      }

      //! @brief copy constructor
      //! @param ASSIGNMENT the assignment to copy
      Assignment( const Assignment< t_Member> &ASSIGNMENT) :
        m_Members( ASSIGNMENT.m_Members)
      {
      }

      //! @brief Clone is the virtual copy constructor
      //! @return pointer to a new Assignment which is a copy of this Assignment
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

      //! @brief GetMembers gives a non-changeable reference to m_Members
      //! @return a non-changeable reference to m_Members
      const util::SiPtrList< const t_Member> &GetMembers() const
      {
        return m_Members;
      }

      //! @brief number of member is assignment
      size_t GetSize() const
      {
        return m_Members.GetSize();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief append a member to the end of the assignment
      //! @param SP_MEMBER an object of t_Member
      void Append( const util::SiPtr< const t_Member> &SP_MEMBER)
      {
        m_Members.Append( SP_MEMBER);
      }

      //! @brief append a list of members to the end of the assignment
      //! @param SP_MEMBER_LIST a list of objects of t_Member
      void Append( const util::SiPtrList< const t_Member> &SP_MEMBER_LIST)
      {
        m_Members.Append( SP_MEMBER_LIST);
      }

      //! @brief compares the SiPtr of the t_Members, returns true if all of them are equal and thus t_Members identical
      //! @param ASSIGNMENT the assignment to compare to *this
      //! @return true if identical
      bool IsIdentical( const Assignment< t_Member> &ASSIGNMENT) const
      {
        return m_Members == ASSIGNMENT.m_Members;
      }

      //! @brief compares t_Members, returns true if all of them are equal (SiPtrs do not need to be identical)
      //! @param ASSIGNMENT the assignment to compare to *this
      //! @return true if t_Members are equal
      bool IsEqual( const Assignment< t_Member> &ASSIGNMENT) const
      {
        // loop over t_Member of *this and ASSIGNMENT in parallel
        for
        (
          typename util::SiPtrList< const t_Member>::const_iterator
            itr( m_Members.Begin()), itr_end( m_Members.End()),
            assignment_itr( ASSIGNMENT.GetMembers().Begin()), assignment_itr_end( ASSIGNMENT.GetMembers().End());
          itr != itr_end && assignment_itr != assignment_itr_end;
          ++itr, ++assignment_itr
        )
        {
          const util::SiPtr< const t_Member> &sp_member( *itr);
          const util::SiPtr< const t_Member> &sp_assignment_member( *assignment_itr);
          // if exactly one ptr is undefined (gap), return false (if both are undef, it's still ok)
          if( !sp_member.IsDefined() || !sp_assignment_member.IsDefined())
          {
            return !sp_member.IsDefined() && !sp_assignment_member.IsDefined();
          }
          if( !( *sp_member == *sp_assignment_member)) // compare t_Member
          {
            return false;
          }
        }

        return true;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief pretty-prints an assignment as string of its member's char IDs (for AAs their one-letter-code)
      //! @return the string
      std::string ToString( const char GAP_CHAR = s_DefaultGapChar) const
      {
        std::string assignment_string;
        for
        (
          typename util::SiPtrList< const t_Member>::const_iterator itr( m_Members.Begin()), itr_end( m_Members.End());
          itr != itr_end;
          ++itr
        )
        {
          char aa( itr->IsDefined() ? GetCharId< t_Member>( **itr) : GAP_CHAR);
          assignment_string.append( 1, aa);
        }
        return assignment_string;
      }

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        BCL_Exit( "not implemented yet", -1);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        OSTREAM << m_Members << '\n';

        return OSTREAM;
      }

    }; // template class Assignment

    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> Assignment< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new Assignment< t_Member>())
    );

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ASSIGNMENT_H_
