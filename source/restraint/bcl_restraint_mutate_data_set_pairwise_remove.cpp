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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_remove.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_data_set_pairwise.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseRemove::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseRemove())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseRemove::GetDefaultScheme()
    {
      static const std::string s_scheme( "remove");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from optional scheme and defaults number of elements to remove to 1
    //! @param SCHEME optional scheme
    MutateDataSetPairwiseRemove::MutateDataSetPairwiseRemove( const std::string &SCHEME) :
      m_SizeRange( 1, 1),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking members
    //! @param MIN the minimum possible number of elements desired to be removed
    //! @param MAX the maximum possible number of elements desired to be removed
    //! @param SCHEME scheme for mutate
    MutateDataSetPairwiseRemove::MutateDataSetPairwiseRemove( const size_t MIN, const size_t MAX, const std::string &SCHEME) :
      m_SizeRange( MIN, MAX),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseRemove
    MutateDataSetPairwiseRemove *MutateDataSetPairwiseRemove::Clone() const
    {
      return new MutateDataSetPairwiseRemove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseRemove::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateDataSetPairwiseRemove::GetScheme() const
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
    //! @param DATA DataSetPairwise of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA_SET
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseRemove::operator()( const DataSetPairwise &DATA) const
    {
      // static empty model
      static util::ShPtr< DataSetPairwise> s_empty_data_set;

      if( DATA.IsEmpty())
      {
        // return skipped move
        return math::MutateResult< DataSetPairwise>( s_empty_data_set, *this);
      }

      util::ShPtr< DataSetPairwise> new_data_set( DATA.Clone());

      // get random number of elements to remove within range
      const size_t num_to_remove( random::GetGlobalRandom().Random( m_SizeRange));

      // remove as many elements as specified
      for( size_t remove( 0); remove < num_to_remove && ( num_to_remove - remove) <= new_data_set->GetSize(); ++remove)
      {
        DataSetPairwise::const_iterator remove_element
        (
          random::GetGlobalRandom().Iterator( new_data_set->Begin(), new_data_set->End(), new_data_set->GetSize())
        );

        new_data_set->Erase( remove_element);
      }

      return math::MutateResult< DataSetPairwise>( new_data_set, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseRemove::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SizeRange, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseRemove::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
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
