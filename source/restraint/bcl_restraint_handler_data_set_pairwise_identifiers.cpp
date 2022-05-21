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
#include "restraint/bcl_restraint_handler_data_set_pairwise_identifiers.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"

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
    const util::SiPtr< const util::ObjectInterface> HandlerDataSetPairwiseIdentifiers::s_Instance
    (
      GetObjectInstances().AddInstance( new HandlerDataSetPairwiseIdentifiers())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    HandlerDataSetPairwiseIdentifiers::HandlerDataSetPairwiseIdentifiers() :
      m_Score( util::GetUndefinedDouble()),
      m_DataSetPairwise()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param SCORE the score of the data set for reading or writing
    HandlerDataSetPairwiseIdentifiers::HandlerDataSetPairwiseIdentifiers( const double SCORE) :
      m_Score( SCORE),
      m_DataSetPairwise()
    {
    }

    //! @brief Clone function
    //! @return pointer to new HandlerDataSetPairwiseIdentifiers
    HandlerDataSetPairwiseIdentifiers *HandlerDataSetPairwiseIdentifiers::Clone() const
    {
      return new HandlerDataSetPairwiseIdentifiers( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerDataSetPairwiseIdentifiers::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the score
    //! @return double which is the score of the dataset
    const double HandlerDataSetPairwiseIdentifiers::GetScore() const
    {
      return m_Score;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads data set from a stream
    //! @param ISTREAM the stream the data set will be read from
    //! @return stream that the data set was read from
    std::istream &HandlerDataSetPairwiseIdentifiers::ReadDataSetPairwise( std::istream &ISTREAM)
    {
      std::string score_id;
      ISTREAM >> score_id >> m_Score;

      // read in the data set
      while( !ISTREAM.eof() && ISTREAM.peek() != std::istream::traits_type::eof())
      {
        DataPairwise data_pair
        (
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( new LocatorCoordinatesFirstSideChainAtom()),
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( new LocatorCoordinatesFirstSideChainAtom())
        );

        // read in current data pair
        data_pair.ReadIdentification( ISTREAM);

        // insert the data pair into the data set
        BCL_MessageDbg( "trying to insert " + data_pair.GetIdentification() + " into " + util::Format()( m_DataSetPairwise));
        BCL_Assert( m_DataSetPairwise.Insert( data_pair).second, "could not insert data pair " + data_pair.GetIdentification());
        BCL_MessageDbg( "successfully inserted " + data_pair.GetIdentification());
        // get next char
        char tmp;
        ISTREAM.get( tmp);
      }

      return ISTREAM;
    }

    //! @brief provides the data set that the handler created
    //! @return const reference to an data set object
    const DataSetPairwise &HandlerDataSetPairwiseIdentifiers::GetDataSetPairwise() const
    {
      return m_DataSetPairwise;
    }

    //! @brief writes restraints to a stream
    //! @param OSTREAM the stream the data set will be written to
    //! @param DATA_SET the data set that will be written to the stream
    //! @return stream the data set were written to
    std::ostream &HandlerDataSetPairwiseIdentifiers::WriteDataSetPairwise
    (
      std::ostream &OSTREAM, const DataSetPairwise &DATA_SET
    ) const
    {
      // print score
      OSTREAM << "score: " << m_Score << '\n';

      // iterate through the sorted data
      for
      (
        DataSetPairwise::const_iterator
        itr( DATA_SET.Begin()), itr_end( DATA_SET.End());
        itr != itr_end; ++itr
      )
      {
        OSTREAM << itr->GetIdentification() << '\n';
      }

      return OSTREAM;
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
    std::istream &HandlerDataSetPairwiseIdentifiers::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Score, ISTREAM);
      io::Serialize::Read( m_DataSetPairwise, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &HandlerDataSetPairwiseIdentifiers::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Score, OSTREAM, INDENT);
      io::Serialize::Write( m_DataSetPairwise, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
