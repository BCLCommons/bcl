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
#include "score/bcl_score_data_set_pairwise_size.h"

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
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseSize::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseSize())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseSize::GetDefaultScheme()
    {
      static const std::string s_scheme( "data_set_size");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseSize::DataSetPairwiseSize( const std::string &SCHEME) :
      m_SizeRange(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking optional scheme
    //! @param SIZE_RANGE_MIN minimum desired size inclusive
    //! @param SIZE_RANGE_MAX maximum desired size inclusive
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseSize::DataSetPairwiseSize
    (
      const size_t &SIZE_RANGE_MIN, const size_t &SIZE_RANGE_MAX, const std::string &SCHEME
    ) :
      m_SizeRange( math::Range< double>( SIZE_RANGE_MIN, SIZE_RANGE_MAX)),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseSize
    DataSetPairwiseSize *DataSetPairwiseSize::Clone() const
    {
      return new DataSetPairwiseSize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseSize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseSize::GetScheme() const
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
    double DataSetPairwiseSize::operator()( const restraint::DataSetPairwise &DATA) const
    {
      const double data_size( DATA.GetSize());

      // true if size is within the desired range
      if( m_SizeRange.IsWithin( data_size))
      {
        return 0;
      }

      // score is the distance outside of the range
      const double score( math::Absolute( data_size - m_SizeRange.GetMiddle()) - ( m_SizeRange.GetWidth() * 0.5));

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseSize::Read( std::istream &ISTREAM)
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
    std::ostream &DataSetPairwiseSize::Write( std::ostream &OSTREAM, const size_t INDENT) const
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

  } // namespace score
  
} // namespace bcl
