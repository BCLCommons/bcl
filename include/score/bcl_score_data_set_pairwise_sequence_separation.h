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

#ifndef BCL_SCORE_DATA_SET_PAIRWISE_SEQUENCE_SEPARATION_H_
#define BCL_SCORE_DATA_SET_PAIRWISE_SEQUENCE_SEPARATION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetPairwiseSequenceSeparation
    //! @brief Implementation of sequence separation term
    //! @details See Kazmier et al. "Algorithm for selection of optimized EPR distance restraints for de novo protein
    //!          structure prediction". Journal of Structural Biology 2011 section 2.1.1.
    //!
    //!          Might need to calculate on per chain basis right now calculates on total size of all chains
    //!
    //! @see @link example_score_data_set_pairwise_sequence_separation.cpp @endlink
    //! @author alexanns
    //! @date May 7, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetPairwiseSequenceSeparation :
      public math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double>
    {

    private:

    //////////
    // data //
    //////////

      //! the sequence size the score will be normalized against
      size_t m_SequenceLength;

      //! the scheme of this mutate
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking optional scheme
      //! @param SCHEME the scheme for this scoring function
      explicit DataSetPairwiseSequenceSeparation( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking scheme
      //! @param SEQUENCE_SIZE the sequence size the score will be normalized against
      //! @param SCHEME the scheme for this scoring function
      explicit DataSetPairwiseSequenceSeparation
      (
        const size_t SEQUENCE_SIZE, const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new DataSetPairwiseSequenceSeparation
      DataSetPairwiseSequenceSeparation *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the score of a data set
      //! @param DATA data set to be scored
      //! @return the score of the current data set
      double operator()( const restraint::DataSetPairwise &DATA) const;

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

      //! @brief determines the sequence separation between two chain and sequence ids
      //! @param CHAIN_ID_A first chain id
      //! @param SEQ_ID_A first seq id
      //! @param CHAIN_ID_B second chain id
      //! @param SEQ_ID_B second seq id
      //! @return size_t which is the sequence separation between the first and second
      static size_t CalculateSequenceSeparation
      (
        const char CHAIN_ID_A, const int SEQ_ID_A, const char CHAIN_ID_B, const int SEQ_ID_B
      );

    }; // class DataSetPairwiseSequenceSeparation

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_DATA_SET_PAIRWISE_SEQUENCE_SEPARATION_H_ 
