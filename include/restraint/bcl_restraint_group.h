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

#ifndef BCL_RESTRAINT_GROUP_H_
#define BCL_RESTRAINT_GROUP_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Group
    //! @brief This is a template class for storing a group of objects
    //!
    //! @see @link example_restraint_group.cpp @endlink
    //! @author alexanns
    //! @date 02/08/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Group :
      public util::SiPtrList< const t_DataType>
    {

    private:

    //////////
    // data //
    //////////

    public:

    /////////////
    // typedef //
    /////////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct a SiPtrList from optional size and single element
      //! @param SIZE optional size, default 0
      //! @param DATA optional single element which is assigned to all stored elements, default is empty SiPtr
      Group< t_DataType>( const size_t SIZE = 0, const util::SiPtr< const t_DataType> &DATA = util::SiPtr< const t_DataType>()) :
        util::SiPtrList< const t_DataType>( SIZE, DATA)
      {
      }

      //! @brief construct from size and pointer to data
      //! @param SIZE number of elements behind the pointer *DATA
      //! @param DATA pointer to the elements
      Group< t_DataType>( const size_t SIZE, t_DataType *const DATA) :
        util::SiPtrList< const t_DataType>( SIZE, DATA)
      {
      }

      //! construct from util::SiPtrList< t_DataType>
      Group< t_DataType>( const util::SiPtrList< const t_DataType> &LIST) :
        util::SiPtrList< const t_DataType>( LIST)
      {
      }

      //! @brief construct from iterator [FIRST, LAST) range
      //! @tparam t_Iterator which is the type of iterator that indicates the range
      //! @param FIRST t_Iterator to the first element to be copied
      //! @param LAST t_Iterator to the first element after the last element to be copied
      template< typename t_Iterator>
      Group< t_DataType>( const t_Iterator &FIRST, const t_Iterator &LAST) :
        util::SiPtrList< const t_DataType>( FIRST, LAST)
      {
      }

      //! @brief construct from two Groups
      //! @param GROUP_A the Group< t_DataType> which will be first in the new Group
      //! @param GROUP_B the Group< t_DataType> which will be second in the new Group
      Group< t_DataType>( const Group< t_DataType> &GROUP_A, const Group< t_DataType> &GROUP_B) :
        util::SiPtrList< const t_DataType>()
      {
        // insert "GROUP_B" at the beginning of the new Group
        util::SiPtrList< const t_DataType>::InsertElements( util::SiPtrList< const t_DataType>::Begin(), GROUP_B);

        // insert "GROUP_A" at the beginning of the new Group, before "GROUP_B"
        util::SiPtrList< const t_DataType>::InsertElements( util::SiPtrList< const t_DataType>::Begin(), GROUP_A);
      }

      //! @brief construct from a Group and a number of empty members
      //! @param GROUP is the Group< t_DataType> which will be inserted in the beginning of the new Group
      //! @param SIZE is the number of empty members the new Group will hold after construction
      Group< t_DataType>( const Group< t_DataType> &GROUP, const size_t SIZE) :
        util::SiPtrList< const t_DataType>( SIZE)
      {
        // insert "GROUP" into the beginning of the new Group before the empty members
        util::SiPtrList< const t_DataType>::InsertElements( util::SiPtrList< const t_DataType>::Begin(), GROUP);
      }

      //! @brief construct from a number of empty members and a Group
      //! @param SIZE is the number of empty members the new Group will hold after construction
      //! @param GROUP is the Group< t_DataType> which will be inserted at the end of the new Group
      Group< t_DataType>( const size_t SIZE, const Group< t_DataType> &GROUP) :
        util::SiPtrList< const t_DataType>( SIZE)
      {
        // insert "GROUP" at the end of the new Group after the empty members
        util::SiPtrList< const t_DataType>::InsertElements( util::SiPtrList< const t_DataType>::End(), GROUP);
      }

      //! @brief copy constructor
      //! @param GROUP the Group to copy
      Group< t_DataType>( const Group< t_DataType> &GROUP) :
        util::SiPtrList< const t_DataType>( GROUP)
      {
      }

      //! virtual copy constructor
      virtual Group< t_DataType> *Clone() const
      {
        return new Group< t_DataType>( *this);
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

    ///////////////
    // operators //
    ///////////////

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
        // write base class
        util::SiPtrList< const t_DataType>::Write( OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

      //! @brief read container from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        // exit since you can't read in a simple pointer
        BCL_Exit( "trying to read in Group, but cannot read in simple pointers", -1);

        // return
        return ISTREAM;
      }

    }; // class Group

  } // namespace restraint
} // namespace bcl

#endif //BCL_RESTRAINT_GROUP_H_
