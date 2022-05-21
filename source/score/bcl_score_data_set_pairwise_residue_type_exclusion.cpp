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
#include "score/bcl_score_data_set_pairwise_residue_type_exclusion.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_data_set_pairwise.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseResidueTypeExclusion::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseResidueTypeExclusion())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseResidueTypeExclusion::GetDefaultScheme()
    {
      static const std::string s_scheme( "aa_type_excl");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DataSetPairwiseResidueTypeExclusion::DataSetPairwiseResidueTypeExclusion( const std::string &SCHEME) :
      m_AATypesToExclude(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking scheme
    //! @param AA_TYPES the types of amino acids which are undesirable
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseResidueTypeExclusion::DataSetPairwiseResidueTypeExclusion
    (
      const storage::Set< biol::AAType> AA_TYPES, const std::string &SCHEME
    ) :
      m_AATypesToExclude( AA_TYPES),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseResidueTypeExclusion
    DataSetPairwiseResidueTypeExclusion *DataSetPairwiseResidueTypeExclusion::Clone() const
    {
      return new DataSetPairwiseResidueTypeExclusion( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseResidueTypeExclusion::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseResidueTypeExclusion::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the score of a data set
    //! @param DATA data set to be scored
    //! @return the score of the current data set
    double DataSetPairwiseResidueTypeExclusion::operator()( const restraint::DataSetPairwise &DATA) const
    {
      double score( 0);

      // iterate through the data set
      for
      (
        restraint::DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // get the aa types of the current residues
        const biol::AAType &aa_type_a( data_itr->First()->GetAAType());
        const biol::AAType &aa_type_b( data_itr->Second()->GetAAType());

        // add the number of times the aa type is found in m_AATypesToExclude to the score
        // this should be either 1 or 0 since the aa types in the set are unique
        score += m_AATypesToExclude.Count( aa_type_a);
        score += m_AATypesToExclude.Count( aa_type_b);
      }

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseResidueTypeExclusion::Read( std::istream &ISTREAM)
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
    std::ostream &DataSetPairwiseResidueTypeExclusion::Write( std::ostream &OSTREAM, const size_t INDENT) const
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

  } // namespace score
  
} // namespace bcl
