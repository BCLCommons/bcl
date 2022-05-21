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

#ifndef BCL_RESTRAINT_ASSIGNMENT_H_
#define BCL_RESTRAINT_ASSIGNMENT_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_group_collection.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Assignment
    //! @brief This is a template class for storing objects which go together/are assigned together for some reason.
    //!
    //! @tparam t_Restraint
    //! @tparam t_GroupIdentifier is the type of object that will be used as the object for creating Groups
    //! @tparam t_GroupMember is the type of object that the Assignment will be storing and putting into Groups
    //! @tparam t_GroupIdentifierCompare is the type of functor that is used to sort the t_GroupIdentifiers as they
    //!           are stored in the Assignment; the default is std::less< t_GroupIdentifier>
    //!
    //! @see @link example_restraint_assignment.cpp @endlink
    //! @author meilerj, alexanns
    //! @date 03/09/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template
    <
      typename t_Restraint, typename t_GroupIdentifier, typename t_GroupMember, typename t_GroupIdentifierCompare
    >
    class Assignment :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      util::ShPtr< t_Restraint> m_Restraint; //!< restraint

      //! "m_GroupCollection" holds the t_GroupMembers and holds them in the appropriate Groups
      GroupCollection< t_GroupIdentifier, t_GroupMember, t_GroupIdentifierCompare> m_GroupCollection;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Assignment() :
        m_Restraint(),
        m_GroupCollection()
      {
      }

      //! @brief construct from a ShPtr to a t_Restraint
      //! @param RESTRAINT is the ShPtr< t_Restraint> object that the Assignment will have
      Assignment( const util::ShPtr< t_Restraint> &RESTRAINT) :
        m_Restraint( RESTRAINT),
        m_GroupCollection()
      {
      }

      //! @brief construct from a GroupCollection
      //! @param GROUP_COLLECTION is the GroupCollection that the Assignment will have
      Assignment( const GroupCollection< t_GroupIdentifier, t_GroupMember, t_GroupIdentifierCompare> &COLLECTION) :
        m_Restraint(),
        m_GroupCollection( COLLECTION)
      {
      }

      //! @brief construct from a ShPtr to a t_Restraint and a GroupCollection
      //! @param RESTRAINT is the ShPtr< t_Restraint> object that the Assignment will have
      //! @param GROUP_COLLECTION is the GroupCollection that the Assignment will have
      Assignment
      (
        const util::ShPtr< t_Restraint>                                                          &RESTRAINT,
        const GroupCollection< t_GroupIdentifier, t_GroupMember, t_GroupIdentifierCompare> &GROUP_COLLECTION
      ) :
        m_Restraint( RESTRAINT),
        m_GroupCollection( GROUP_COLLECTION)
      {
      }

      //! @brief Clone is the virtual copy constructor
      //! @return returns a pointer to a new Assignment which is a copy of this Assignment
      virtual Assignment< t_Restraint, t_GroupIdentifier, t_GroupMember, t_GroupIdentifierCompare> *Clone() const
      {
        return new Assignment( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief GetRestraint gives a non-changeable reference to "m_Restraint"
      //! @return return a non-changeable reference to the ShPtr< t_Restraint> which is "m_Restraint"
      virtual util::ShPtr< t_Restraint> const &GetRestraint() const
      {
        return m_Restraint;
      }

      //! @brief SetRestraint changes "m_Restraint" to a new ShPtr< t_Restraint>
      //! @param RESTRAINT is the ShPtr< t_Restraint> "m_Restraint" will be changed to
      virtual void SetRestraint( const util::ShPtr< t_Restraint> &RESTRAINT)
      {
        m_Restraint = RESTRAINT;
      }

      //! @brief GetGroupCollection gives a non-changeable reference to "m_GroupCollection"
      //! @return return a non-changeable reference to the GroupCollection which is "m_GroupCollection"
      virtual GroupCollection< t_GroupIdentifier, t_GroupMember, t_GroupIdentifierCompare> const &GetGroupCollection() const
      {
        return m_GroupCollection;
      }

      //! @brief SetGroupCollection changes "m_GroupCollection" to a new GroupCollection
      //! @param GROUP_COLLECTION is the GroupCollection "m_GroupCollection" will be changed to
      virtual void
      SetGroupCollection
      (
        const GroupCollection< t_GroupIdentifier, t_GroupMember, t_GroupIdentifierCompare> &GROUP_COLLECTION
      )
      {
        m_GroupCollection = GROUP_COLLECTION;
      }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! write assignment to std::ostream
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        OSTREAM << m_Restraint << '\n';
        OSTREAM << m_GroupCollection << '\n';

        return OSTREAM;
      }

      //! read assignment from std::istream
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        BCL_Exit( "not implemented yet", -1);
        return ISTREAM;
      }

    }; // class Assignment

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ASSIGNMENT_H_
