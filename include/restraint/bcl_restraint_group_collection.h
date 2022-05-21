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

#ifndef BCL_RESTRAINT_GROUP_COLLECTION_H_
#define BCL_RESTRAINT_GROUP_COLLECTION_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_group.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GroupCollection
    //! @brief This is a template class for storing a collection of Groups that are separated by identifiers.
    //!
    //! @tparam t_GroupIdentifier is the type of object that will be used to identify and denote the Groups
    //! @tparam t_GroupMember is the type of objects that will be held in the Groups
    //! @tparam t_GroupIdentifierCompare is the type of functor that will be used to sort the GroupCollection by
    //!           the t_GroupIdentifiers: the default is the std::less_than
    //!
    //! @see @link example_restraint_group_collection.cpp @endlink
    //! @author alexanns
    //! @date 02/09/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_GroupIdentifier, typename t_GroupMember, typename t_GroupIdentifierCompare>
    class GroupCollection :
      public storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>
    {

    private:

    //////////
    // data //
    //////////

    public:

    /////////////
    // typedef //
    /////////////

      //! iterator
      typedef typename storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::iterator iterator;

      //! const iterator
      typedef typename storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::const_iterator const_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      GroupCollection() :
        storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>()
      {
      }

      //! @brief construct from an initial group identifier and Group
      //! @param GROUP_IDENTIFIER is the t_GroupIdentifier which is associated with the Group to be added
      //! @param GROUP is the Group< t_GroupMember> which will be added to the new GroupCollection
      GroupCollection
      (
        const t_GroupIdentifier     &GROUP_IDENTIFIER,
        const Group< t_GroupMember> &GROUP
      ) :
        storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>()
      {
        // insert the "GROUP_IDENTIFIER" and "GROUP" into the new GroupCollection
        storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::Insert( std::pair< t_GroupIdentifier, Group< t_GroupMember> >( GROUP_IDENTIFIER, GROUP));
      }

      //! virtual copy constructor
      virtual GroupCollection< t_GroupIdentifier, t_GroupMember, t_GroupIdentifierCompare> *Clone() const
      {
        return new GroupCollection< t_GroupIdentifier, t_GroupMember, t_GroupIdentifierCompare>( *this);
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

    ////////////////
    // operations //
    ////////////////

      using storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::Insert;

      //! @brief Insert takes a t_GroupIdentifier and a Group and inserts it into the GroupCollection
      //! @param GROUP_IDENTIFIER is the t_GroupIdentifier which is associated with the Group to be added
      //! @param GROUP is the Group< t_GroupMember> which will be added to the new GroupCollection
      //! @return a std pair with an iterator pointing to the place of insertion and a bool indicating success or not
      //! if the element already existed then the iterator points to the previous existance and the bool is false
      std::pair< iterator, bool>
      Insert( const t_GroupIdentifier &GROUP_IDENTIFIER, const Group< t_GroupMember> &GROUP)
      {
        // insert the "GROUP_IDENTIFIER" and "GROUP" into the GroupCollection and return result
        return storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::Insert( std::pair< t_GroupIdentifier, Group< t_GroupMember> >( GROUP_IDENTIFIER, GROUP));
      }

      //! @brief SingleGroupDepth gives the size of the Group indicated by an iterator
      //! @param ITERATOR is the iterator which indicates the Group whose size is of interest
      //! @return size_t which indicates the number of members in the Group indicated by "ITERATOR"
      virtual size_t
      SingleGroupDepth
      (
        const const_iterator ITERATOR
      ) const
      {
        return ITERATOR->second.GetSize();
      }

      //! @brief TotalDepth gives the total number of members within all Groups of the Assignment
      //! @return size_t which is the sum of the size of each Group within the Assignment
      virtual size_t TotalDepth() const
      {
        // create "total_depth" to hold the total number of members within the Group Collection
        size_t total_depth( 0);

        // iterate over the GroupCollection
        for
        (
          typename GroupCollection< t_GroupIdentifier, t_GroupMember, t_GroupIdentifierCompare>::const_iterator
            group_collection_iterator( storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::Begin()),
            group_collection_iterator_end( storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::End());
          group_collection_iterator != group_collection_iterator_end;
          ++group_collection_iterator
        )
        {
          // add the depth of the Group currently denoted by "group_collection_iterator" to "total_depth"
          total_depth += SingleGroupDepth( group_collection_iterator);
        }

          // return "total_depth"
          return total_depth;
      }

      //! @brief GetOverallNthGroupMember gives an iterator to the member of the Assignment specified by its
      //!        sequential numerical position. The range of possible arguments is [0, TotalDepth()).
      //! @param NTH_GROUP_MEMBER is the size_t which denotes the member of interest
      //! @return returns a const_iterator which denotes the member of interest
      virtual typename Group< t_GroupMember>::const_iterator
      GetOverallNthGroupMember( size_t NTH_GROUP_MEMBER) const
      {
        // make sure that the Nth group member exists in the Assignment
        BCL_Assert
        (
          NTH_GROUP_MEMBER < TotalDepth(),
          "trying to get the " + util::Format()( NTH_GROUP_MEMBER) + " th group member but only " +
          util::Format()( TotalDepth()) + " members exist"
         );

        // create const_iterator "group_collection_iterator" and initialize to the beginning of this GroupCollection
        const_iterator group_collection_iterator( storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::Begin());

        size_t current_group_depth( SingleGroupDepth( group_collection_iterator));

        // create size_t "number_members" which will keep the running sum of members found in each Group
        size_t number_members( current_group_depth);

        // create const_iterator on a Group "group_iterator"
        typename Group< t_GroupMember>::const_iterator group_iterator;

        // if "number_members" is greater than "NTH_GROUP_MEMBER" then the Group denoted by "group_collection_iterator"
        // (the first Group of the GroupCollection) contains the member of interest
        if( number_members > NTH_GROUP_MEMBER)
        {
          // set "group_iterator" to the beginning of the Group denoted by "group_collection_iterator"
          group_iterator = group_collection_iterator->second.Begin();

          // move "group_iterator" to the correct position in the Group
          for( size_t x( 0); x < NTH_GROUP_MEMBER; ++x)
          {
            ++group_iterator;
          }

          // return "group_iterator"
          return group_iterator;
        }

        // create size_t "previous_sum" which will hold the sum of members up to before the Group denoted by
        // "group_collection_iterator"
        size_t previous_sum( 0);

        // iterate through the GroupCollection and sum up "number_members" until it is greater than "NTH_GROUP_MEMBER"
        while( number_members <= NTH_GROUP_MEMBER)
        {
          // move "group_collection_iterator" to the next Group
          ++group_collection_iterator;

          previous_sum += current_group_depth;

          current_group_depth = SingleGroupDepth( group_collection_iterator);

          // check if the Group currently indicated by "group_collection_iterator" has the desired member
          number_members += current_group_depth;
        }

        // set "group_iterator" to the beginning of the Group denoted by "group_collection_iterator"
        group_iterator = group_collection_iterator->second.Begin();

        // move "group_iterator" to the correct position in the Group
        for( size_t x( 0); x < ( NTH_GROUP_MEMBER - previous_sum); ++x)
        {
          ++group_iterator;
        }

        // return "group_iterator"
        return group_iterator;
      }

      //! @brief CollectGroupMembers puts all of the t_GroupMembers of this Assignment into a single Group
      //! @return gives a Group< t_GroupMember> that contains all the members of all Groups from the GroupCollection
      virtual Group< t_GroupMember> CollectAllGroupMembers() const
      {
        // create Group which will hold all of the t_GroupMembes of all Groups in the Assignment
        Group< t_GroupMember> big_group;

        // iterate through the GroupCollection and add the t_GroupMembers of all Groups to "big_group"
        for
        (
          typename GroupCollection< t_GroupIdentifier, t_GroupMember, t_GroupIdentifierCompare>::const_iterator
            group_collection_iterator( storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::Begin()),
            group_collection_iterator_end( storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::End());
          group_collection_iterator != group_collection_iterator_end;
          ++group_collection_iterator
        )
        {
          // add the group denoted by "group_collection_iterator" to the end of "big_group"
          big_group.Append( group_collection_iterator->second);
        }

          // return "big_group" which as all the t_GroupMembers of this Assignment's GroupCollections
          return big_group;
      }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief write container to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write m_Groups
        storage::Map< t_GroupIdentifier, Group< t_GroupMember>, t_GroupIdentifierCompare>::Write( OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

      //! @brief read container from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        // exit since you can't read in a simple pointers of Groups
        BCL_Exit( "trying to read in GroupCollection, but cannot read in simple pointers of Groups", -1);

        // return
        return ISTREAM;
      }

    }; // class GroupCollection

  } // namespace restraint
} // namespace bcl

#endif //BCL_RESTRAINT_GROUP_COLLECTION_H_
