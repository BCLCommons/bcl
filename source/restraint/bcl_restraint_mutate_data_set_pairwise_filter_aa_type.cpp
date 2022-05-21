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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_aa_type.h"

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
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterAAType::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterAAType())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterAAType::GetDefaultScheme()
    {
      static const std::string s_scheme( "filter_aa_type_excl");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateDataSetPairwiseFilterAAType::MutateDataSetPairwiseFilterAAType( const std::string &SCHEME) :
      m_AATypesToExclude(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking scheme
    //! @param AA_TYPES the types of amino acids which are undesirable
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterAAType::MutateDataSetPairwiseFilterAAType
    (
      const storage::Set< biol::AAType> AA_TYPES, const std::string &SCHEME
    ) :
      m_AATypesToExclude( AA_TYPES),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterAAType
    MutateDataSetPairwiseFilterAAType *MutateDataSetPairwiseFilterAAType::Clone() const
    {
      return new MutateDataSetPairwiseFilterAAType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterAAType::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterAAType::GetScheme() const
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
    //! @param SEQUENCE DataSetPairwise of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA_SET
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseFilterAAType::operator()( const DataSetPairwise &SEQUENCE) const
    {
      static util::ShPtr< DataSetPairwise> s_empty;

      if( SEQUENCE.IsEmpty() || m_AATypesToExclude.IsEmpty())
      {
        BCL_MessageStd
        (
          "data set or AATypes to exclude is empty. data set size : " +
          util::Format()( SEQUENCE.GetSize()) + " AATypes to exclude size is : " +
          util::Format()( m_AATypesToExclude.GetSize())
        );
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // the dataset filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      // iterate through the data set to filter out data pairs that don't meet exposure criteria
      for
      (
        DataSetPairwise::const_iterator data_itr( SEQUENCE.Begin()), data_itr_end( SEQUENCE.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // get the aa type of the locators
        const biol::AAType &aa_type_a( data_itr->First()->GetAAType());
        const biol::AAType &aa_type_b( data_itr->Second()->GetAAType());

        // true if neither of the aa types are in m_AATypesToExclude
        if( !m_AATypesToExclude.Contains( aa_type_a) && !m_AATypesToExclude.Contains( aa_type_b))
        {
          // add data pair to filtered_data
          filtered_data->Insert( *data_itr);
        }
      }

      BCL_MessageStd
      (
        "MutateDataSetPairwiseFilterAAType num unique locators " +
        util::Format()( filtered_data->GetUniqueDataPoints().GetSize())
      );

      BCL_MessageStd( "end MutateDataSetPairwiseFilterAAType::operator()");
      return math::MutateResult< DataSetPairwise>( filtered_data, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseFilterAAType::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AATypesToExclude, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseFilterAAType::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AATypesToExclude, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
