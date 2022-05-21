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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_add.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseAdd::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseAdd())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseAdd::GetDefaultScheme()
    {
      static const std::string s_scheme( "add");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from optional scheme
    //! @param SCHEME optional scheme
    MutateDataSetPairwiseAdd::MutateDataSetPairwiseAdd( const std::string &SCHEME) :
      m_CompleteDataSet(),
      m_SizeRange(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking sequence
    //! @param POOL_DATA_SET pool of possible data points to add
    //! @param SCHEME optional scheme
    MutateDataSetPairwiseAdd::MutateDataSetPairwiseAdd
    (
      const util::ShPtr< DataSetPairwise> &POOL_DATA_SET, const std::string &SCHEME
    ) :
      m_CompleteDataSet( POOL_DATA_SET),
      m_SizeRange( 1, 1),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking sequence
    //! @param POOL_DATA_SET pool of possible data points to add
    //! @param MIN min possible number that will be added
    //! @param MAX max possible number that will be added
    //! @param SCHEME optional scheme
    MutateDataSetPairwiseAdd::MutateDataSetPairwiseAdd
    (
      const util::ShPtr< DataSetPairwise> &POOL_DATA_SET, const size_t MIN, const size_t MAX, const std::string &SCHEME
    ) :
      m_CompleteDataSet( POOL_DATA_SET),
      m_SizeRange( MIN, MAX),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseAdd
    MutateDataSetPairwiseAdd *MutateDataSetPairwiseAdd::Clone() const
    {
      return new MutateDataSetPairwiseAdd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseAdd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateDataSetPairwiseAdd::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an DataSetPairwise and returning a mutated object
    //! @param DATA_SET DataSetPairwise of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA_SET
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseAdd::operator()( const DataSetPairwise &DATA_SET) const
    {
//      // will hold the data points that are available to add to to the current dataset
//      std::vector< DataPairwise> difference;
//
//      BCL_Assert( m_CompleteDataSet.IsDefined(), "complete data set pointer is not defined");
//
//      std::set_difference
//      (
//        m_CompleteDataSet->Begin(), m_CompleteDataSet->End(),
//        DATA_SET.Begin(), DATA_SET.End(),
//        std::inserter( difference, difference.begin())
//      );

      // static empty data set
      static util::ShPtr< DataSetPairwise> s_empty_data_set;

      // true if no available data pairs to add
      if( m_CompleteDataSet->IsEmpty()) //difference.empty())
      {
        BCL_MessageDbg( "difference is empty");
        // return skipped move
        return math::MutateResult< DataSetPairwise>( s_empty_data_set, *this);
      }

      util::ShPtr< DataSetPairwise> new_data_set( DATA_SET.Clone());

//      // randomly shuffle the difference data pairs
//      std::random_shuffle( difference.begin(), difference.end());

      // get random number of elements to add within range but no more than are in difference
      const size_t num_to_add( std::min( random::GetGlobalRandom().Random( m_SizeRange), m_CompleteDataSet->GetSize())); // difference.size()));

      // add as many elements are necessary to the new data set
      for( size_t num_added( 0); num_added < num_to_add; ++num_added)
      {
//        // add the as many elements as needed
//        new_data_set->InsertElements( difference.begin(), difference.begin() + num_to_add);

        // get random iterator to one of the data pairs in the pool of possibles
        DataSetPairwise::const_iterator random_itr
        (
          random::GetGlobalRandom().Iterator
          (
            m_CompleteDataSet->Begin(), m_CompleteDataSet->End(), m_CompleteDataSet->GetSize()
          )
        );
        BCL_Assert( random_itr != m_CompleteDataSet->End(), "random iterator at end of total dataset");

        // insert random data pair from the complete data set into the new data set
        new_data_set->Insert( *random_itr);
      }
      BCL_MessageDbg
      (
        "MutateDataSetPairwiseAdd data set size is " + util::Format()( new_data_set->GetSize())
      );
      BCL_Assert( new_data_set.IsDefined(), "new_data_set pointer is not defined");
      return math::MutateResult< DataSetPairwise>( new_data_set, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseAdd::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_CompleteDataSet, ISTREAM);
      io::Serialize::Read( m_SizeRange, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseAdd::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_CompleteDataSet, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SizeRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
