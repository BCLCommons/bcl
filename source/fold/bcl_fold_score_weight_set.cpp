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
#include "fold/bcl_fold_score_weight_set.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ScoreWeightSet::s_Instance
    (
      GetObjectInstances().AddInstance( new ScoreWeightSet())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ScoreWeightSet::ScoreWeightSet() :
        m_WeightMap()
    {
    }

    //! @brief constructor from a table
    //! @param TABLE that contains the weights
    ScoreWeightSet::ScoreWeightSet( const storage::Table< double> &TABLE) :
        m_WeightMap()
    {
      InitializeFromTable( TABLE);
    }

    //! @brief Clone function
    //! @return pointer to new ScoreWeightSet
    ScoreWeightSet *ScoreWeightSet::Clone() const
    {
      return new ScoreWeightSet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ScoreWeightSet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ScoreWeightSet::GetAlias() const
    {
      static const std::string s_name( "ScoreWeights");
      return s_name;
    }

    //! @brief returns the weight map
    //! @return weight map
    const storage::Map< Score, double> &ScoreWeightSet::GetWeightMap() const
    {
      return m_WeightMap;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get weight for a given score
    //! @param SCORE Score of interest
    //! @return weight for the given score
    double ScoreWeightSet::GetWeight( const Score &SCORE)
    {
      return m_WeightMap[ SCORE];
    }

    //! @brief sets the weight of the given score
    //! @param SCORE enum of the score
    //! @param WEIGHT Weight to be assigned to score
    //! @return true is score is defined and weight was set; false otherwise
    bool ScoreWeightSet::SetWeight
    (
      const Score &SCORE,
      const double WEIGHT
    )
    {
      // make sure the a score is defined
      if( !SCORE.IsDefined())
      {
        BCL_MessageVrb( "undefined score supplied!");
        return false;
      }

      // update the weight
      m_WeightMap[ SCORE] = WEIGHT;

      // end
      return true;
    }

    //! @brief sets the weight of the given score
    //! @param SCORE_NAME name of the score
    //! @param WEIGHT Weight to be assigned to score
    //! @return true is score name is a valid score name and weight was set; false otherwise
    bool ScoreWeightSet::SetWeight
    (
      const std::string &SCORE_NAME,
      const double WEIGHT
    )
    {
      // construct the score enum from the name and try to set the weight
      return SetWeight( Score( SCORE_NAME), WEIGHT);
    }

    //! @brief resets all weights to zero
    void ScoreWeightSet::Reset()
    {
      // iterate over the map
      for
      (
        storage::Map< Score, double>::iterator map_itr( m_WeightMap.Begin()), map_itr_end( m_WeightMap.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // set the weight to zero
        map_itr->second = 0.0;
      }
    }

    //! @brief constructs the scores and returns it
    //! @return constructed score
    util::ShPtr< score::ProteinModelScoreSum> ScoreWeightSet::ConstructScoreSum() const
    {
      // construct the score
      util::ShPtr< score::ProteinModelScoreSum> scoresum( new score::ProteinModelScoreSum( m_WeightMap));

//      // true if an ensemble of models is being folded
//      if( fold::ProtocolEnsemble::GetFlagEnsembleSize()->GetFlag())
//      {
//        BCL_MessageCrt( "setting score sum to use ensemble score sum");
//
//        // set scoresum pointer to a protein ensemble score sum
//        scoresum = util::ShPtr< score::ProteinModelScoreSum>( new score::ProteinEnsembleScoreSum( m_WeightMap));
//      }

      // return the scoresum object
      return scoresum;
    }

    //! @brief creates a table from scoring weight set
    //! @return table with the scoring weight set
    storage::Table< double> ScoreWeightSet::CreateTable() const
    {
      // create a vector to hold column names
      storage::Vector< std::string> column_names;
      // create a vector to hold the weightst
      storage::Vector< double> weights;

      // iterate over the map
      for
      (
        storage::Map< Score, double>::const_iterator map_itr( m_WeightMap.Begin()), map_itr_end( m_WeightMap.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // insert the name and the weight
        column_names.PushBack( map_itr->first.GetName());
        weights.PushBack( map_itr->second);
      }

      // create a table
      storage::Table< double> table( column_names);
      // insert a row
      table.InsertRow( "weights", weights);

      // end
      return table;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ScoreWeightSet::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_WeightMap, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ScoreWeightSet::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_WeightMap, OSTREAM);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initialize this ScoreWeightSet from a table
    //! @param TABLE Table that contains the weightset
    void ScoreWeightSet::InitializeFromTable( const storage::Table< double> &TABLE)
    {
      // if the table does not have row weights
      if( !TABLE.HasRow( "weights"))
      {
        BCL_Assert
        (
          TABLE.GetHeader().HasColumn( "weights"),
          "The given score weightset table has no row or column named \"weights\""
        );

        // transpose table
        storage::Table< double> transposed_table( TABLE.GetTransposedTable());

        // initialize with the transposed table
        return InitializeFromTable( transposed_table);
      }

      // get the map from the row
      storage::Map< std::string, double> string_map( TABLE[ "weights"].ConvertToMap());

      // iterate over the map
      for
      (
        storage::Map< std::string, double>::const_iterator map_itr( string_map.Begin()), map_itr_end( string_map.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // get the score enum with the corresponding name
        const Score &score( GetScores().GetEnumFromName( map_itr->first));

        // make sure there is a score with the given string
        BCL_Assert
        (
          score != GetScores().e_Undefined,
          "There is no score with the name \"" + map_itr->first + "\" among:\n" +
          util::Format()( storage::Vector< std::string>( GetScores().GetEnumStrings()))
        );

        // set the weight for this score in the weight map
        m_WeightMap[ score] = map_itr->second;
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ScoreWeightSet::GetSerializer() const
    {
      io::Serializer params;
      params.SetClassDescription( "scores and weights for folding");
      params.AddInitializer
      (
        "",
        "the map with scores and weights",
        io::Serialization::GetAgent( &m_WeightMap)
      );
      return params;
    }

  } // namespace fold
} // namespace bcl
