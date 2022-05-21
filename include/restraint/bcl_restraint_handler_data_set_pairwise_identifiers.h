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

#ifndef BCL_RESTRAINT_HANDLER_DATA_SET_PAIRWISE_IDENTIFIERS_H_
#define BCL_RESTRAINT_HANDLER_DATA_SET_PAIRWISE_IDENTIFIERS_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_data_set_pairwise.h"
#include "bcl_restraint_handler_data_set_pairwise_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerDataSetPairwiseIdentifiers
    //! @brief allows reading and writing data set in specific format
    //! @details the format is
    //!
    //! score: 3
    //! 'A'  32  "UsePDBID"  0  "<=>"  'A'  48  "UsePDBID"  0
    //!
    //! @see @link example_restraint_handler_data_set_pairwise_identifiers.cpp @endlink
    //! @author alexanns
    //! @date Oct 28, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API HandlerDataSetPairwiseIdentifiers :
      public HandlerDataSetPairwiseInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the score of the data set for reading or writing
      double m_Score;

      //! the data set this handler has read in
      DataSetPairwise m_DataSetPairwise;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      HandlerDataSetPairwiseIdentifiers();

      //! @brief constructor taking member variable parameters
      //! @param SCORE the score of the data set for reading or writing
      HandlerDataSetPairwiseIdentifiers( const double SCORE);

      //! @brief Clone function
      //! @return pointer to new HandlerDataSetPairwiseIdentifiers
      HandlerDataSetPairwiseIdentifiers *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the score
      //! @return double which is the score of the dataset
      const double GetScore() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief reads data set from a stream
      //! @param ISTREAM the stream the data set will be read from
      //! @return stream that the data set was read from
      std::istream &ReadDataSetPairwise( std::istream &ISTREAM);

      //! @brief provides the data set that the handler created
      //! @return const reference to an data set object
      const DataSetPairwise &GetDataSetPairwise() const;

      //! @brief writes restraints to a stream
      //! @param OSTREAM the stream the data set will be written to
      //! @param DATA_SET the data set that will be written to the stream
      //! @return stream the data set were written to
      std::ostream &WriteDataSetPairwise( std::ostream &OSTREAM, const DataSetPairwise &DATA_SET) const;

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

    private:

    }; // class HandlerDataSetPairwiseIdentifiers

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_HANDLER_DATA_SET_PAIRWISE_IDENTIFIERS_H_ 
