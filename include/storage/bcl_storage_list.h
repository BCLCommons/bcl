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

#ifndef BCL_STORAGE_LIST_H_
#define BCL_STORAGE_LIST_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialize.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically
#include <list>

namespace bcl
{
  namespace storage
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class List
    //! @brief This is a template class for List Containers.
    //! @details Attributes: (copied from http://www.sgi.com/tech/stl/List.html)
    //! - double linked list (allows forward and backwards traversal)
    //! - constant time insertion and removal of elements at any position
    //! - iterators are invalidated only if the iterator was pointing to an element that was removed
    //! - no random access to elements
    //!
    //! @tparam t_DataType indicates the type of data that will be held
    //!
    //! @see @link example_storage_list.cpp @endlink
    //! @author meilerj, heinzes1
    //! @date 08/03/2004, updated 08/30/2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class List :
      public util::ObjectInterface
    {

    private:

      std::list< t_DataType> m_Data;

    //////////
    // data //
    //////////

    public:

      //! typedef for iterator
      typedef typename std::list< t_DataType>::iterator               iterator;
      //! typedef for const_iterator
      typedef typename std::list< t_DataType>::const_iterator         const_iterator;
      //! typedef for reverse_iterator
      typedef typename std::list< t_DataType>::reverse_iterator       reverse_iterator;
      //! typedef for const_reverse_iterator
      typedef typename std::list< t_DataType>::const_reverse_iterator const_reverse_iterator;
      //! typedef for const_reference
      typedef typename std::list< t_DataType>::reference              reference;
      //! typedef for const_reference
      typedef typename std::list< t_DataType>::const_reference        const_reference;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      List() :
        m_Data()
      {
      }

      //! @brief construct a List from optional size and single element
      //! @param SIZE optional size, default 0
      //! @param VALUE optional single element which is assigned to all stored elements, default is t_DataType()
      List( const size_t SIZE, const t_DataType &VALUE = t_DataType()) :
        m_Data( SIZE, VALUE)
      {
      }

      //! @brief construct a List from number of size and pointer to data
      //! @param SIZE number of elements behind the pointer *DATA
      //! @param DATA pointer to the elements
      List( const size_t SIZE, const t_DataType *DATA) :
        m_Data( SIZE)
      {
        for( iterator itr( m_Data.begin()), itr_end( m_Data.end()); itr != itr_end; ++itr, ++DATA)
        {
          *itr = *( DATA);
        }
      }

      //! @brief construct a List from another List
      //! @param LIST the List to copy
      List( const List< t_DataType> &LIST) :
        m_Data( LIST.m_Data)
      {
      }

      //! @brief move a List from another List
      //! @param LIST the List to copy
      List( List< t_DataType> && LIST) :
        m_Data( std::move( LIST.m_Data))
      {
      }

      //! @brief construct a List from iterator [FIRST, LAST) range
      //! create a List from the given range [BEGIN, END) excluding END
      //! @tparam t_Iterator which is the type of iterator that indicates the range
      //! @param BEGIN t_Iterator to the first element to be copied
      //! @param END t_Iterator to the first element after the last copied element
      template< typename t_Iterator>
      List( const t_Iterator &BEGIN, const t_Iterator &END) :
        m_Data( BEGIN, END)
      {
      }

      //! @brief virtual copy constructor
      //! @return pointer to a copy of the actual object
      List< t_DataType> *Clone() const
      {
        return new List< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief GetClassIdentifier returns class name of the object
      //! @return returns std::string which is the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns size of the container
      //! @return size, i.e. number of elements stored
      size_t GetSize() const
      {
        return m_Data.size();
      }

      //! @brief return the maximal size of the container
      //! @return maximum size, i.e. maximal number of elements to store
      size_t MaxSize() const
      {
        return m_Data.max_size();
      }

      //! @brief returns a const reference to the used stl container
      //! @return const reference to the internal stl container
      std::list< t_DataType> const &InternalData() const
      {
        return m_Data;
      }

      //! @brief returns a changeable reference to the used stl container
      //! @return changeable reference to the internal stl container
      std::list< t_DataType> &InternalData()
      {
        return m_Data;
      }

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      iterator Begin()
      {
        return m_Data.begin();
      }

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      const_iterator Begin() const
      {
        return m_Data.begin();
      }

      //! @brief return changeable iterator to last element
      //! @return changeable iterator to last element
      iterator Last()
      {
        BCL_Assert( !IsEmpty(), "Cannot return iterator to last element from empty sequence container");
        return --m_Data.end();
      }

      //! @brief return const iterator to last element
      //! @return const iterator to last element
      const_iterator Last() const
      {
        BCL_Assert( !IsEmpty(), "Cannot return iterator to last element from empty sequence container");
        return --m_Data.end();
      }

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      iterator End()
      {
        return m_Data.end();
      }

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      const_iterator End() const
      {
        return m_Data.end();
      }

      //! @brief return iterator to reverse begin
      //! @return reverse_iterator pointing to the beginning of the reversed container, i.e. the last element
      reverse_iterator ReverseBegin()
      {
        return m_Data.rbegin();
      }

      //! @brief return const_iterator to reverse begin
      //! @return const_reverse_iterator pointing to the beginning of the reversed container
      const_reverse_iterator ReverseBegin() const
      {
        return m_Data.rbegin();
      }

      //! @brief return iterator to reverse end
      //! @return reverse_iterator pointing to the end of the reversed container, i.e. behind the first element
      reverse_iterator ReverseEnd()
      {
        return m_Data.rend();
      }

      //! @brief return const_iterator to reverse end
      //! @return const_reverse_iterator pointing to the end of the reversed container
      const_reverse_iterator ReverseEnd() const
      {
        return m_Data.rend();
      }

      //! @brief return const reference to first element
      //! @return const reference to first element of type t_DataType
      const t_DataType &FirstElement() const
      {
        return m_Data.front();
      }

      //! @brief return a changeable reference to first element
      t_DataType &FirstElement()
      {
        return m_Data.front();
      }

      //! @brief return const reference to last element
      //! @return const reference to last element of type t_DataType
      const t_DataType &LastElement() const
      {
        return m_Data.back();
      }

      //! @brief return changeable reference to last element
      //! @return changeable reference to last element of type t_DataType
      t_DataType &LastElement()
      {
        return m_Data.back();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief checks whether container is empty
      //! @return if the container is empty or not
      bool IsEmpty() const
      {
        return m_Data.empty();
      }

      //! @brief delete all elements
      void Reset()
      {
        m_Data.clear();
      }

      //! @brief insert ELEMENT into the container
      //! @param ELEMENT an object of t_DataType that is inserted
      void InsertElement( const t_DataType &ELEMENT)
      {
        m_Data.push_back( ELEMENT);
      }

      //! @brief insert ELEMENT before ITR into the container
      //! @param ITR iterator on the Sequence Container where ELEMENT is inserted
      //! @param ELEMENT an object of t_DataType that is inserted
      template< typename t_InputIterator>
      void InsertElement( const t_InputIterator ITR, const t_DataType &ELEMENT)
      {
        m_Data.insert( ITR, ELEMENT);
      }

      //! @brief insert all elements of argument CONTAINER
      //! @param ITR_POS iterator
      //! @param BEGIN start location
      //! @param END end location
      template< typename t_InputIterator>
      void InsertElements( iterator ITR_POS, t_InputIterator BEGIN, t_InputIterator END)
      {
        m_Data.insert( ITR_POS, BEGIN, END);
      }

      //! @brief insert all elements of argument CONTAINER
      //! @param ITR_POS iterator
      //! @param LIST list
      void InsertElements( iterator ITR_POS, const List< t_DataType> &LIST)
      {
        m_Data.insert( ITR_POS, LIST.Begin(), LIST.End());
      }

      //! @brief delete single element at ITR
      //! @return iterator pointing to the element immediuately following the deleted one
      //! @param ITR iterator pointing to the element that will be destroyed
      iterator Remove( iterator ITR)
      {
        return ITR == End() ? End() : m_Data.erase( ITR);
      }

      //! @brief delete a range of elements [FIRST, LAST)
      //! @param ITR_FIRST first element of range to be deleted
      //! @param ITR_LAST last element, in the range, will be kept
      void Remove( iterator ITR_FIRST, iterator ITR_LAST)
      {
        m_Data.erase( ITR_FIRST, ITR_LAST);
      }

      //! @brief delete single element at ITR
      //! @param ITR iterator pointing to the element that will be destroyed
      void RemoveElement( iterator ITR)
      {
        ITR == End() ? End() : m_Data.erase( ITR);
      }

      //! @brief append a ELEMENT to the end of *this
      //! @param ELEMENT an object of t_DataType that is appended
      void Append( const t_DataType &ELEMENT)
      {
        m_Data.push_back( ELEMENT);
      }

      //! @brief append a ELEMENT to the end of *this
      //! @param ELEMENT an object of t_DataType that is appended
      void Append( t_DataType && ELEMENT)
      {
        m_Data.push_back( std::move( ELEMENT));
      }

      //! @brief append a LIST to the end of *this
      //! @param LIST list of t_DataType that is appended
      void Append( const List< t_DataType> &LIST)
      {
        m_Data.insert( End(), LIST.Begin(), LIST.End());
      }

      //! @brief insert all elements of argument CONTAINER
      //! @param BEGIN start location
      //! @param END end location
      template< typename t_InputIterator>
      void Append( t_InputIterator BEGIN, t_InputIterator END)
      {
        m_Data.insert( End(), BEGIN, END);
      }

      //! @brief append a ELEMENT to the end of *this
      //! @param ELEMENT an object of t_DataType that is appended
      void Prepend( const t_DataType &ELEMENT)
      {
        m_Data.push_front( ELEMENT);
      }

      //! @brief append a LIST to the end of *this
      //! @param LIST list of t_DataType that is appended
      void Prepend( const List< t_DataType> &LIST)
      {
        m_Data.insert( Begin(), LIST.Begin(), LIST.End());
      }

      //! @brief insert all elements of argument CONTAINER
      //! @param BEGIN first element
      //! @param END last element
      template< typename t_InputIterator>
      void Prepend( t_InputIterator BEGIN, t_InputIterator END)
      {
        m_Data.insert( Begin(), BEGIN, END);
      }

      //! @brief inserts ELEMENT to the end
      //! @param ELEMENT element to be inserted, default element is t_DataType()
      void PushBack( const t_DataType &ELEMENT = t_DataType())
      {
        m_Data.push_back( ELEMENT);
      }

      //! @brief inserts ELEMENT to the end
      //! @param ELEMENT element to be inserted, default element is t_DataType()
      void PushBack( t_DataType && ELEMENT)
      {
        m_Data.push_back( std::move( ELEMENT));
      }

      //! @brief removes last element
      void PopBack()
      {
        m_Data.pop_back();
      }

      //! @brief resizes the sequence to the length NUMBER_ELEMENTS, eventually fill with ELEMENTs
      //! @param NUMBER_ELEMENTS number of elements to resize to
      //! @param ELEMENT optional element of t_DataType to fill the Sequence Container with
      void Resize( const size_t NUMBER_ELEMENTS, const t_DataType &ELEMENT = t_DataType())
      {
        m_Data.resize( NUMBER_ELEMENTS, ELEMENT);
      }

      //! @brief set all elements of the sequence to one given value X
      //! @param ELEMENT the element of type t_DataType which will be assigned to all elements, default is t_DataType()
      void SetAllElements( const t_DataType &ELEMENT = t_DataType())
      {
        for( iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          ( *itr) = ELEMENT;
        }
      }

      //! @brief PushFront inserts ELEMENT to the front
      //! @param ELEMENT t_DataType element to be inserted, default element is t_DataType()
      void PushFront( const t_DataType &ELEMENT = t_DataType())
      {
        m_Data.push_front( ELEMENT);
      }

      //! @brief PopFront removes first element
      void PopFront()
      {
        m_Data.pop_front();
      }

      //! @brief all elements of List X are inserted in *this before position ITR
      //! and removed from X; *this must not be X
      //!
      //! This function splices the List X before position ITR in the List *this.
      //! X is empty after this operation.
      //! The other splice() functions are NOT implemented because of compile errors:
      //!   template< typename t_Iterator> void Splice( t_Iterator ITR, List< t_DataType> &X, t_Iterator X_ITR)
      //! and
      //!   template< typename t_Iterator> void
      //!   Splice( t_Iterator ITR, List< t_DataType> &X, t_Iterator X_FIRST, t_Iterator X_LAST)
      //! Using the std::list functions leads to functions which either do
      //! nothing or in case of an invalid X_ITR in an endless loop.
      //!
      //! @tparam t_Iterator type of iterator which will be used to denote where to start splicing
      //! @param ITR t_Iterator denoting position before the elements of X are inserted
      //! @param X List< t_DataType> with elements for insertion in *this
      template< typename t_Iterator>
      void Splice( t_Iterator ITR, List< t_DataType> &X)
      {
        m_Data.splice( ITR, X.InternalData());
      }

      //! @brief Sort sorts the list according to a binary predicate
      //! in approx N log N, where N is the length of the list
      //! @tparam t_BinaryPredicate type of binary predicate that will be used for comparison
      //! @param COMPARABLE t_BinaryPredicate used to compare two elements of the list
      template< typename t_BinaryPredicate>
      void Sort( t_BinaryPredicate COMPARABLE)
      {
        m_Data.sort( COMPARABLE);
      }

      //! @brief Reverse reverses the order of elements in the list
      void Reverse()
      {
        m_Data.reverse();
      }

      //! @brief Unique removes all but the first element in every consecutive group of equal elements
      //! @tparam t_BinaryPredicate type of binary predicate that will be used for comparison
      //! @param COMPARABLE t_BinaryPredicate binary predicate to compare two elements of the list
      template< typename t_BinaryPredicate>
      void Unique( t_BinaryPredicate COMPARABLE)
      {
        m_Data.unique( COMPARABLE);
      }

      //! @brief Merge merges two sorted Lists (*this and another) together according to binary predicate COMPARABLE;
      //! The second List is empty after the merge
      //! @tparam t_BinaryPredicate type of binary predicate that will be used for comparison
      //! @param X List< t_DataType> is the sorted List to be merged with *this
      //! @param COMPARABLE t_BinaryPredicate used to compare two elements of the list
      template< typename t_BinaryPredicate>
      void Merge( List< t_DataType> &X, t_BinaryPredicate COMPARABLE)
      {
        m_Data.merge( X.InternalData(), COMPARABLE);
      }

      //! @brief Extracts a number of elements as a block starting from a given position
      //! @param POS size_t to determine the position of the element to be extracted
      //! @param NUMBER number of elements to be extracted
      List< t_DataType> ExtractElements( size_t POS, size_t NUMBER = 1)
      {
        // initialize two lists to store the remaining data and the extracted portion
        List< t_DataType> extract, remain;
        size_t actual_position( 0);

        // iterate over the list
        for
        (
          const_iterator
            itr( Begin()), itr_end( End());
          itr != itr_end; ++itr, ++actual_position
        )
        {
          // if before or after the requested fraction skip and put in the remaining list
          if( actual_position < POS || actual_position >= POS + NUMBER)
          {
            remain.PushBack( *itr);
          }
          // otherwise add to the extracted list
          else
          {
            extract.PushBack( *itr);
          }
        }

        // assign the remaining list to this object ( "this" was remain + extract)
        ( *this) = remain;

        // return the extracted
        return extract;
      }

      //! @brief Extracts one element from a given position
      //! @param POS size_t to determine the position of the element to be extracted
      t_DataType ExtractElement( size_t POS)
      {
        return ExtractElements( POS, 1).FirstElement();
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief equal operator
      List< t_DataType> &operator =( const List< t_DataType> &LIST)
      {
        if( this != &LIST)
        {
          m_Data = LIST.m_Data;
        }
        return *this;
      }

      //! @brief move assignment operator
      List< t_DataType> &operator =( List< t_DataType> && LIST)
      {
        if( this != &LIST)
        {
          m_Data = std::move( LIST.m_Data);
        }
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write container to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write size of container
        io::Serialize::Write( GetSize(), OSTREAM, INDENT) << '\n';

        // write each element of container to new line
        for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          io::Serialize::Write( *itr, OSTREAM, INDENT) << '\n';
        }

        // end
        return OSTREAM;
      }

      //! @brief read container from io::IFStream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // create int "size" for holding the size of the container
        size_t size;

        // read in size of container
        io::Serialize::Read( size, ISTREAM);

        // match container size to size of incoming data
        Resize( size);

        // read in each element from stream
        for( iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          BCL_Assert( io::Serialize::Read( *itr, ISTREAM), "Error reading element!");
        }

        // return
        return ISTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief create from two elements
      static List< t_DataType> Create( const t_DataType &DATA_A, const t_DataType &DATA_B)
      {
        List< t_DataType> list;
        list.PushBack( DATA_A);
        list.PushBack( DATA_B);
        return list;
      }

    }; //end template class List

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> List< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new List< t_DataType>())
    );

    //! @brief operator == checks if two containers are the same
    //! @param LIST_A first container
    //! @param LIST_B second container
    //! @return returns boolean value
    template< typename t_DataType>
    inline bool operator ==
    (
      const List< t_DataType> &LIST_A,
      const List< t_DataType> &LIST_B
    )
    {
      return LIST_A.InternalData() == LIST_B.InternalData();
    }

  } // namespace storage
} // namespace bcl

#endif //BCL_STORAGE_LIST_H_
