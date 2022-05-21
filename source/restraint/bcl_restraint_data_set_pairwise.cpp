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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "restraint/bcl_restraint_data_set_pairwise.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwise::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwise())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DataSetPairwise::DataSetPairwise() :
      m_DataSet()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Sequence
    DataSetPairwise *DataSetPairwise::Clone() const
    {
      return new DataSetPairwise( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwise::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return iterator on begin
    //! @return iterator pointing to the beginning of the container, i.e. the first element
    DataSetPairwise::iterator DataSetPairwise::Begin()
    {
      return m_DataSet.Begin();
    }

    //! @brief return const_iterator on begin
    //! @return const_iterator pointing to the beginning of the container, i.e. the first element
    DataSetPairwise::const_iterator DataSetPairwise::Begin() const
    {
      return m_DataSet.Begin();
    }

    //! @brief return iterator on end
    //! @return iterator pointing to the end of the container, i.e. behind the last element
    DataSetPairwise::iterator DataSetPairwise::End()
    {
      return m_DataSet.End();
    }

    //! @brief return const_iterator on end
    //! @return const_iterator pointing to the end of the container, i.e. behind the last element
    DataSetPairwise::const_iterator DataSetPairwise::End() const
    {
      return m_DataSet.End();
    }

    //! @brief return iterator to reverse begin
    //! @return reverse_iterator pointing to the beginning of the reversed container, i.e. the last element
    DataSetPairwise::reverse_iterator DataSetPairwise::ReverseBegin()
    {
      return m_DataSet.ReverseBegin();
    }

    //! @brief return const_iterator to reverse begin
    //! @return const_reverse_iterator pointing to the beginning of the reversed container
    DataSetPairwise::const_reverse_iterator DataSetPairwise::ReverseBegin() const
    {
      return m_DataSet.ReverseBegin();
    }

    //! @brief return iterator to reverse end
    //! @return reverse_iterator pointing to the end of the reversed container, i.e. behind the first element
    DataSetPairwise::reverse_iterator DataSetPairwise::ReverseEnd()
    {
      return m_DataSet.ReverseEnd();
    }

    //! @brief return const_iterator to reverse end
    //! @return const_reverse_iterator pointing to the end of the reversed container
    DataSetPairwise::const_reverse_iterator DataSetPairwise::ReverseEnd() const
    {
      return m_DataSet.ReverseEnd();
    }

    //! @brief indicates if there is nothing in the data set or not
    //! @return bool true if dataset is empty -false otherwise
    bool DataSetPairwise::IsEmpty() const
    {
      return m_DataSet.IsEmpty();
    }

    //! @brief gives the number of elements in the data set
    //! @return size_t indicating the number of elements in the data set
    size_t DataSetPairwise::GetSize() const
    {
      return m_DataSet.GetSize();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief AddPair function for adding a pair of atoms to the dataset
    //! @param LOCATOR_A first of the two atoms to be added to the dataset
    //! @param LOCATOR_B second of the atoms to be paired with the first
    //! @return pair of bool true if insertion was successful - false otherwise and iterator pointing to insertion place
    std::pair< DataSetPairwise::iterator, bool> DataSetPairwise::Insert
    (
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_A,
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_B
    )
    {
      // make an empty pairwise data
      DataPairwise data;

      // set the data
      data.Set( LOCATOR_A, LOCATOR_B);

      // try to insert the data and return the result
      return std::pair< iterator, bool>( m_DataSet.Insert( data));
    }

    //! @brief AddPair function for adding a pair of atoms to the dataset
    //! @param DATA data pair to be added
    //! @return pair of bool  true if insertion was successful - false otherwise and iterator pointing to place of insertion
    std::pair< DataSetPairwise::iterator, bool> DataSetPairwise::Insert
    (
      const DataPairwise &DATA
    )
    {
      return std::pair< iterator, bool>( m_DataSet.Insert( DATA));
    }

    //! @brief Erase function for removing a data pair from the data set
    //! @param DATA_PAIR the data points to be removed
    //! @return bool true if DATA_PAIR was removed - false otherwise
    bool DataSetPairwise::Erase( const DataPairwise &DATA_PAIR)
    {
      return m_DataSet.Erase( DATA_PAIR);
    }

    //! @brief gives all of the single points that are unique in the dataset
    //! @return set which has the individual data points that are unique from the dataset
    storage::Set
    <
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
      assemble::LocatorAtomCoordinatesInterface::PtrLessThan
    > DataSetPairwise::GetUniqueDataPoints() const
    {
      // get the list of all data points
      storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > data_points;

      // iterate through the data set to fill up all the data points
      for
      (
        DataSetPairwise::const_iterator data_itr( Begin()), data_itr_end( End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // try to insert the data points, only unique points will be able to be inserted
        data_points.Insert( data_itr->First());
        data_points.Insert( data_itr->Second());
      }

      return data_points;
    }

    //! @brief gives all of the single points in the dataset - leaves in duplicates
    //! @return set which has the individual data points leaving in duplicates
    util::ShPtrList< assemble::LocatorAtomCoordinatesInterface> DataSetPairwise::GetDataPoints() const
    {
      util::ShPtrList< assemble::LocatorAtomCoordinatesInterface> data_points;

      // iterate through the data set to fill up all the data points
      for
      (
        DataSetPairwise::const_iterator data_itr( Begin()), data_itr_end( End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // try to insert the data points
        data_points.PushBack( data_itr->First());
        data_points.PushBack( data_itr->Second());
      }

      // sort the data points by sequence
      data_points.Sort( assemble::LocatorAtomCoordinatesInterface::PtrLessThan());

      return data_points;
    }

    //! @brief Find takes a key and returns an iterator to the element denoted by key, or an off the end iterator
    //! @param KEY the key for which an iterator is desired
    //! @return returns an iterator to the element
    DataSetPairwise::iterator DataSetPairwise::Find( const DataPairwise &KEY)
    {
      return m_DataSet.Find( KEY);
    }

    //! @brief Find takes a key and returns a const_iterator to the element denoted by key, or off the end iterator
    //! @param KEY the key for which a const_iterator is desired
    //! @return returns an iterator to the element
    DataSetPairwise::const_iterator DataSetPairwise::Find( const DataPairwise &KEY) const
    {
      return m_DataSet.Find( KEY);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwise::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DataSet, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwise::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DataSet, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief creates a complete pairwise dataset from vector of aabases
    //! @param SEQUENCE the sequence for which a complete data set will be created
    //! @return DataSetPairwise which has the complete dataset for SEQUENCE
    util::ShPtr< DataSetPairwise> DataSetPairwise::GetCompleteDataSet( const util::ShPtrVector< biol::AABase> &SEQUENCE)
    {
      util::ShPtr< DataSetPairwise> dataset( new DataSetPairwise());

      // iterate through the sequence
      for
      (
        biol::AASequence::const_iterator aa_itr_a( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr_a != aa_itr_end;
        ++aa_itr_a
      )
      {
        // make locator a
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_a
        (
          assemble::LocatorAtomCoordinatesInterface::CreateLocator< LocatorCoordinatesFirstSideChainAtom>
          (
            **aa_itr_a
          )
        );

        // iterate through the sequence again
        for
        (
          biol::AASequence::const_iterator aa_itr_b( ++biol::AASequence::const_iterator( aa_itr_a));
          aa_itr_b != aa_itr_end;
          ++aa_itr_b
        )
        {
          // make locator b
          const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_b
          (
            assemble::LocatorAtomCoordinatesInterface::CreateLocator< LocatorCoordinatesFirstSideChainAtom>
            (
              **aa_itr_b
            )
          );

          // try to insert locator a and b into the dataset
          const bool inserted( dataset->Insert( locator_a, locator_b).second);

          // make sure the pair could be inserted
          BCL_Assert
          (
            inserted, "could not insert " + locator_a->GetIdentification() + " and " + locator_b->GetIdentification()
          );
        }
      }

      BCL_MessageDbg( "complete dataset size is " + util::Format()( dataset->GetSize()));

      // return the data set
      return dataset;
    }

    //! @brief gets the intersection out of a list of data sets
    //! @param DATA_SETS the data sets whose intersection will be determined
    //! @return a data set that is the intersection of all that have been passed
    DataSetPairwise DataSetPairwise::GetIntersection( const storage::List< DataSetPairwise> &DATA_SETS)
    {
      BCL_Assert( !DATA_SETS.IsEmpty(), "DATA_SETS is empty");

      // get first data set
      DataSetPairwise start_data_set( DATA_SETS.FirstElement());

      // will hold the current intersection
      std::set< DataPairwise> intersection;

      // iterate through the data sets
      for
      (
        storage::List< DataSetPairwise>::const_iterator
          set_itr( ++DATA_SETS.Begin()), set_itr_end( DATA_SETS.End());
        set_itr != set_itr_end; ++set_itr
      )
      {
        // the current set
        const DataSetPairwise &current_set( *set_itr);

        // get the current intersection
        std::set_intersection
        (
          start_data_set.Begin(),      start_data_set.End(),
          current_set.Begin(),         current_set.End(),
          std::inserter( intersection, intersection.begin())
        );

        BCL_MessageDbg( "intersection size is " + util::Format()( intersection.size()));

        // set to the current intersection
        start_data_set = DataSetPairwise( intersection.begin(), intersection.end());

        // clear the current intersection
        intersection.clear();
      }

      // return the intersection of all data sets
      return start_data_set;
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
    )
    {
      // get the list of mutates
      const util::SiPtrList< const math::MutateInterface< DataSetPairwise> > &mutates( MUTATE_RESULT.GetNodes());

      util::SiPtrList< const math::MutateInterface< DataSetPairwise> > all_mutates;

      util::ShPtr< DataSetPairwise> mutated_dataset( MUTATE_RESULT.GetArgument());
      BCL_Assert( mutated_dataset.IsDefined(), "mutated_dataset pointer is not defined");

      // repeat the list of mutates SCALAR number of times
      for( double repeat( 0); repeat < SCALAR - 1; ++repeat)
      {
        // iterate through the mutates
        for
        (
          util::SiPtrList< const math::MutateInterface< DataSetPairwise> >::const_iterator
            mutate_itr( mutates.Begin()), mutate_itr_end( mutates.End());
          mutate_itr != mutate_itr_end;
          ++mutate_itr
        )
        {
          // apply the current mutate
          const math::MutateResult< DataSetPairwise> current_mutate_result( ( *mutate_itr)->operator()( *mutated_dataset));

          mutated_dataset = current_mutate_result.GetArgument();
          BCL_Assert( mutated_dataset.IsDefined(), "mutated_dataset pointer is not defined");

          // add mutates to all_mutates
          all_mutates.Prepend( current_mutate_result.GetNodes());
        }
      }

      // make the final mutate result with final mutated dataset and all the mutates that led to it
      const math::MutateResult< DataSetPairwise> mutate_result( mutated_dataset, all_mutates);
      BCL_Assert( mutated_dataset.IsDefined(), "mutated_dataset pointer is not defined");

      // add the
      return mutate_result;
    }

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
    )
    {
      return operator *( SCALAR, MUTATE_RESULT);
    }

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
    )
    {
      BCL_MessageStd( "plus");
      util::ShPtr< DataSetPairwise> intersection;
      if( LHS.GetArgument().IsDefined() && RHS.GetArgument().IsDefined())
      {
        BCL_MessageStd( "DataSetPairwise::GetIntersection");
        storage::List< DataSetPairwise> results( 1, *LHS.GetArgument());
        results.PushBack( *RHS.GetArgument());
        intersection = util::ShPtr< DataSetPairwise>( DataSetPairwise::GetIntersection( results).Clone());
      }
      else if( LHS.GetArgument().IsDefined())
      {
        BCL_MessageStd( "LHS.GetArgument().IsDefined()");
        intersection = LHS.GetArgument();
      }
      else if( RHS.GetArgument().IsDefined())
      {
        BCL_MessageStd( "RHS.GetArgument().IsDefined()");
        intersection = RHS.GetArgument();
      }
      util::SiPtrList< const math::MutateInterface< DataSetPairwise> > nodes( LHS.GetNodes());
      nodes.Prepend( RHS.GetNodes());
      BCL_Assert( intersection.IsDefined(), "data set pointer is not defined");
      const math::MutateResult< DataSetPairwise> result( intersection, nodes);

      return result;
    }

  } // namespace restraint

} // namespace bcl
