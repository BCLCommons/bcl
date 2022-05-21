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

#ifndef BCL_RESTRAINT_DATA_SET_PAIRWISE_H_
#define BCL_RESTRAINT_DATA_SET_PAIRWISE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_data_pairwise.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetPairwise
    //! @brief This class is for representing a set of data points that are made up of pairs
    //! @details  multiple data pairs
    //!
    //! @see @link example_restraint_data_set_pairwise.cpp @endlink
    //! @author alexanns
    //! @date May 6, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetPairwise :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the set of data that is present
      storage::Set< DataPairwise> m_DataSet;

    public:

      //! typedef for iterator
      typedef storage::Set< DataPairwise>::iterator               iterator;

      //! typedef for const_iterator
      typedef storage::Set< DataPairwise>::const_iterator         const_iterator;

      //! typedef for reverse_iterator
      typedef storage::Set< DataPairwise>::reverse_iterator       reverse_iterator;

      //! typedef for const_reverse_iterator
      typedef storage::Set< DataPairwise>::const_reverse_iterator const_reverse_iterator;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DataSetPairwise();

      //! @brief construct a DataSetPairwise from iterator [FIRST, LAST) range
      //! @tparam t_Iterator which is the type of iterator that indicates the range
      //! @param FIRST t_Iterator to the first element to be copied
      //! @param LAST t_Iterator to the first element after the last copied element
      template< typename t_Iterator>
      DataSetPairwise( t_Iterator FIRST, t_Iterator LAST) :
        m_DataSet( FIRST, LAST)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Sequence
      DataSetPairwise *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      iterator Begin();

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      const_iterator Begin() const;

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      iterator End();

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      const_iterator End() const;

      //! @brief return iterator to reverse begin
      //! @return reverse_iterator pointing to the beginning of the reversed container, i.e. the last element
      reverse_iterator ReverseBegin();

      //! @brief return const_iterator to reverse begin
      //! @return const_reverse_iterator pointing to the beginning of the reversed container
      const_reverse_iterator ReverseBegin() const;

      //! @brief return iterator to reverse end
      //! @return reverse_iterator pointing to the end of the reversed container, i.e. behind the first element
      reverse_iterator ReverseEnd();

      //! @brief return const_iterator to reverse end
      //! @return const_reverse_iterator pointing to the end of the reversed container
      const_reverse_iterator ReverseEnd() const;

      //! @brief indicates if there is nothing in the data set or not
      //! @return bool true if dataset is empty -false otherwise
      bool IsEmpty() const;

      //! @brief gives the number of elements in the data set
      //! @return size_t indicating the number of elements in the data set
      size_t GetSize() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief AddPair function for adding a pair of atoms to the dataset
      //! @param LOCATOR_A first of the two atoms to be added to the dataset
      //! @param LOCATOR_B second of the atoms to be paired with the first
      //! @return pair of bool  true if insertion was successful - false otherwise and iterator pointing to place of insertion
      std::pair< iterator, bool> Insert
      (
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_A,
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_B
      );

      //! @brief AddPair function for adding a pair of atoms to the dataset
      //! @param DATA data pair to be added
      //! @return pair of bool  true if insertion was successful - false otherwise and iterator pointing to place of insertion
      std::pair< iterator, bool> Insert
      (
        const DataPairwise &DATA
      );

      //! @brief InsertElements inserts a number of elements into the container [INCLUSIVE, NOT_INCLUSIVE)
      //! @param ITR_A iterator denoting the first element to be inserted
      //! @param ITR_B iterator denoting the end of the elements to be inserted (not included in the insertions)
      template< typename t_Iterator>
      void InsertElements( t_Iterator ITR_A, t_Iterator ITR_B)
      {
        m_DataSet.InsertElements( ITR_A, ITR_B);
      }

      //! @brief Erase function for removing a data pair from the data set
      //! @param DATA_PAIR the data points to be removed
      //! @return bool true if DATA_PAIR was removed - false otherwise
      bool Erase( const DataPairwise &DATA_PAIR);

      //! @brief Erase function for removing a data pair from the data set at position given by iterator
      //! @param ITR the iterator type indicating the position to be removed
      template< typename t_Iterator>
      void Erase( t_Iterator ITR);

      //! @brief gives all of the single points that are unique in the dataset
      //! @return set which has the individual data points that are unique from the dataset
      storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > GetUniqueDataPoints() const;

      //! @brief gives all of the single points in the dataset - leaves in duplicates
      //! @return set which has the individual data points leaving in duplicates
      util::ShPtrList< assemble::LocatorAtomCoordinatesInterface> GetDataPoints() const;

      //! @brief Find takes a key and returns an iterator to the element denoted by key, or an off the end iterator
      //! @param KEY the key for which an iterator is desired
      //! @return returns an iterator to the element
      iterator Find( const DataPairwise &KEY);

      //! @brief Find takes a key and returns a const_iterator to the element denoted by key, or off the end iterator
      //! @param KEY the key for which a const_iterator is desired
      //! @return returns an iterator to the element
      const_iterator Find( const DataPairwise &KEY) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief creates a complete pairwise dataset from vector of aabases
      //! @param SEQUENCE the sequence for which a complete data set will be created
      //! @return DataSetPairwise which has the complete dataset for SEQUENCE
      static util::ShPtr< DataSetPairwise> GetCompleteDataSet( const util::ShPtrVector< biol::AABase> &SEQUENCE);

      //! @brief gets the intersection out of a list of data sets
      //! @param DATA_SETS the data sets whose intersection will be determined
      //! @return a data set that is the intersection of all that have been passed
      static DataSetPairwise GetIntersection( const storage::List< DataSetPairwise> &DATA_SETS);

    }; // class DataSetPairwise

  ////////////////
  // operations //
  ////////////////

    //! @brief Erase function for removing a data pair from the data set at position given by iterator
    //! @param ITR the iterator type indicating the position to be removed
    template< typename t_Iterator>
    void DataSetPairwise::Erase( t_Iterator ITR)
    {
      return m_DataSet.RemoveElement( ITR);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator * to multiply mutate result of dataset pairwise with a scalar
    //!        repeats the list of mutates nodes on the argument the number of times indicated by the scalar
    //! @param SCALAR scalar to multiply mutate result by
    //! @param MUTATE_RESULT result that will be multiplied
    //! @return new mutate result that has been multiplied by the scalar
    math::MutateResult< DataSetPairwise>
    operator *
    (
      const double &SCALAR,
      const math::MutateResult< DataSetPairwise> &MUTATE_RESULT
    );

    //! @brief operator * to multiply mutate result of dataset pairwise with a scalar
    //!        repeats the list of mutates nodes on the argument the number of times indicated by the scalar
    //! @param SCALAR scalar to multiply mutate result by
    //! @param MUTATE_RESULT result that will be multiplied
    //! @return new mutate result that has been multiplied by the scalar
    math::MutateResult< DataSetPairwise>
    operator *
    (
      const math::MutateResult< DataSetPairwise> &MUTATE_RESULT,
      const double &SCALAR
    );

    //! @brief operator + to add two mutate results of dataset pairwise
    //!        takes the intersection of the two mutate results and combines the mutates
    //! @param LHS first argument
    //! @param RHS second argument
    //! @return new mutate result that is the sum of LHS and RHS
    math::MutateResult< DataSetPairwise>
    operator +
    (
      const math::MutateResult< DataSetPairwise> &LHS,
      const math::MutateResult< DataSetPairwise> &RHS
    );

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_DATA_SET_PAIRWISE_H_ 
