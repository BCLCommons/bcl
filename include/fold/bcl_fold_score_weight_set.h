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

#ifndef BCL_FOLD_SCORE_WEIGHT_SET_H_
#define BCL_FOLD_SCORE_WEIGHT_SET_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_fold_scores.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScoreWeightSet
    //! @brief Convenience class that stores the weightset for scoring function components and constructs the score sum
    //! @details This class stores a map that associates each Score enum with a weight. Through this class you can
    //! change the weights for classes, add new weights and construct the final scoring function ready to be plugged
    //! into an optimizer.
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Nov 19, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreWeightSet :
      public util::SerializableInterface
    {

    private:

    //////////
    // data //
    //////////

      //! map of scores and associated weights
      storage::Map< Score, double> m_WeightMap;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ScoreWeightSet();

      //! @brief constructor from a table
      //! @param TABLE that contains the weights
      ScoreWeightSet( const storage::Table< double> &TABLE);

      //! @brief Clone function
      //! @return pointer to new ScoreWeightSet
      ScoreWeightSet *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief returns the weight map
      //! @return weight map
      const storage::Map< Score, double> &GetWeightMap() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get weight for a given score
      //! @param SCORE Score of interest
      //! @return weight for the given score
      double GetWeight( const Score &SCORE);

      //! @brief sets the weight of the given score
      //! @param SCORE enum of the score
      //! @param WEIGHT Weight to be assigned to score
      //! @return true is score is defined and weight was set; false otherwise
      bool SetWeight
      (
        const Score &SCORE,
        const double WEIGHT
      );

      //! @brief sets the weight of the given score
      //! @param SCORE_NAME name of the score
      //! @param WEIGHT Weight to be asssigned to score
      //! @return true is score name is a valid score name and weight was set; false otherwise
      bool SetWeight
      (
        const std::string &SCORE_NAME,
        const double WEIGHT
      );

      //! @brief constructs the scores and returns it
      //! @return constructed score
      util::ShPtr< score::ProteinModelScoreSum> ConstructScoreSum() const;

      //! @brief creates a table from scoring weight set
      //! @return table with the scoring weight set
      storage::Table< double> CreateTable() const;

      //! @brief resets all weights to zero
      void Reset();

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

      //! @brief initialize this ScoreWeightSet from a table
      //! @param TABLE Table that contains the weightset
      void InitializeFromTable( const storage::Table< double> &TABLE);

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ScoreWeightSet

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_SCORE_WEIGHT_SET_H_
