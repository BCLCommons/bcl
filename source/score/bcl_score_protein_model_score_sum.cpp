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
#include "score/bcl_score_protein_model_score_sum.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_multiplier.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ProteinModelScoreSum::s_Instance
    (
      // while this class, in principle, could be added to the util::Enumerated< ProteinModel> set,
      // math::SumFunctionMixin already satisfies all functionality covered by score::ProteinModel, so there
      // would be no benefit. This class just extends those methods with alternative formatting for non-interface
      // functions
      GetObjectInstances().AddInstance( new ProteinModelScoreSum())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelScoreSum::ProteinModelScoreSum( const std::string &SCHEME) :
      math::SumFunctionMixin< ProteinModel>( SCHEME)
    {
    }

    //! @brief constructor from a map of functions and weights
    //! @param SCORE_WEIGHT_MAP map of scores and corresponding weights
    ProteinModelScoreSum::ProteinModelScoreSum
    (
      const storage::Map< util::ShPtr< ProteinModel>, double> &SCORE_WEIGHT_MAP
    ) :
      math::SumFunctionMixin< ProteinModel>()
    {
      AddScoresWithWeights( SCORE_WEIGHT_MAP);
    }

    //! @brief constructor from a map of Score enums and weights
    //! @param SCORE_WEIGHT_MAP map of scores and corresponding weights
    ProteinModelScoreSum::ProteinModelScoreSum
    (
      const storage::Map< fold::Score, double> &SCORE_WEIGHT_MAP
    ) :
      math::SumFunctionMixin< ProteinModel>()
    {
      // iterate over weights map provided
      for
      (
        storage::Map< fold::Score, double>::const_iterator
          map_itr( SCORE_WEIGHT_MAP.Begin()), map_itr_end( SCORE_WEIGHT_MAP.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // if the weight is equal to 0
        if( map_itr->second == double( 0.0))
        {
          // skip this one
          BCL_MessageStd
          (
            "Weight is equal to 0, therefore not adding the following score " + map_itr->first.GetName()
          );
        }
        // otherwise
        else
        {
          // add the score
          NewOperand( **map_itr->first, map_itr->second);
        }
      }
    }

    //! @brief constructor from a map of functions and from a weight set
    //! @param SCORE_MAP map of scoring functions to be used
    //! @param WEIGHT_SET map of function schemes and corresponding weights
    ProteinModelScoreSum::ProteinModelScoreSum
    (
      const storage::Map< std::string, util::ShPtr< ProteinModel> > &SCORE_MAP,
      const storage::Map< std::string, double> &WEIGHT_SET
    ) :
      math::SumFunctionMixin< ProteinModel>()
    {
      AddScoresWithWeights( SCORE_MAP, WEIGHT_SET);
    }

    //! @brief virtual copy constructor
    ProteinModelScoreSum *ProteinModelScoreSum::Clone() const
    {
      return new ProteinModelScoreSum( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelScoreSum::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns a vector of string that has readable schemes of the individual functions
    //! @return the vector of strings that has readable schemes of the individual functions
    storage::Vector< std::string> ProteinModelScoreSum::GetReadableFunctionSchemes() const
    {
      // get the combined score map
      const storage::Map< ProteinModel::TypeEnum, storage::Map< std::string, double> > score_map( CombineSimilarScores());

      // add some padding to create subheadings
      const std::string padding( "  ");

      // initialize vector to store the function names
      storage::Vector< std::string> function_names;

      // get the structure based scores
      const storage::Map< ProteinModel::TypeEnum, storage::Map< std::string, double> >::const_iterator struct_find
      (
        score_map.Find( ProteinModel::e_Structure)
      );
      if( struct_find != score_map.End())
      {
        function_names.PushBack( ProteinModel::GetTypeName( ProteinModel::e_Structure));
        for
        (
          storage::Map< std::string, double>::const_iterator score_itr( struct_find->second.Begin()),
            score_itr_end( struct_find->second.End());
          score_itr != score_itr_end; ++score_itr
        )
        {
          function_names.PushBack( padding + score_itr->first);
        }
      }

      // get the sequence based scores
      const storage::Map< ProteinModel::TypeEnum, storage::Map< std::string, double> >::const_iterator seq_find
      (
        score_map.Find( ProteinModel::e_Sequence)
      );
      if( struct_find != score_map.End())
      {
        function_names.PushBack( ProteinModel::GetTypeName( ProteinModel::e_Sequence));
        for
        (
          storage::Map< std::string, double>::const_iterator score_itr( seq_find->second.Begin()),
            score_itr_end( seq_find->second.End());
          score_itr != score_itr_end; ++score_itr
        )
        {
          function_names.PushBack( padding + score_itr->first);
        }
      }

      // push sum column
      function_names.PushBack( "Score (total)");

      // end
      return function_names;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &ProteinModelScoreSum::GetAlias() const
    {
      static const std::string s_name( "ProteinModelScoreSum");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief creates a table with individual function values, weights and weighted values
    //! the table has scores as the columns
    //! @param PROTEIN_MODEL ProteinModel to be used for calculating the function value
    //! @return a table with individual functions and their weighted sums
    storage::Table< double> ProteinModelScoreSum::CreateValueTableHorizontal
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // cast a pointer to the multiplier data
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
      );

      // set boolean to see if this is multimer
      assemble::ProteinModel multimer_model;
      const bool has_multimers( sp_multiplier.IsDefined());
      if( has_multimers)
      {
        multimer_model = sp_multiplier->operator ()( PROTEIN_MODEL);
      }

      // initialize vectors to store values and allocate memory
      storage::Vector< double> weights_vector, values_vector, weighted_values_vector;
      weights_vector.AllocateMemory( m_Functions.GetSize() + 2);
      values_vector.AllocateMemory( m_Functions.GetSize() + 2);
      weighted_values_vector.AllocateMemory( m_Functions.GetSize() + 2);

      // initialize the sums
      double sum_weights( 0.0), sum_values( 0.0), sum_weighted_values( 0.0);

      // initialize vector to store the function names
      storage::Vector< std::string> function_names;

      // iterate to collect the scores
      for
      (
        const_iterator function_itr( m_Functions.Begin()), function_itr_end( m_Functions.End());
        function_itr != function_itr_end;
        ++function_itr
      )
      {
        // insert the name for this function
        function_names.PushBack( function_itr->Second()->GetScheme());

        // calculate the value for the argument using this function
        const double this_weight( function_itr->First());
        double this_value( 0.0);
        if( has_multimers)
        {
          this_value = ( function_itr->Second()->operator()( multimer_model));
        }
        else
        {
          this_value = ( function_itr->Second()->operator()( PROTEIN_MODEL));
        }
        const double this_weighted_value( this_value * function_itr->First());

        // record the value, weight and the weighted value
        weights_vector.PushBack( this_weight);
        values_vector.PushBack( this_value);
        weighted_values_vector.PushBack( this_weighted_value);

        // update the sums
        sum_weights += this_weight;
        sum_values += this_value;
        sum_weighted_values += this_weighted_value;
      }

      // add the sum column
      function_names.PushBack( "sum");

      // add the summed values
      weights_vector.PushBack( sum_weights);
      values_vector.PushBack( sum_values);
      weighted_values_vector.PushBack( sum_weighted_values);

      // Insert the rows
      util::ShPtr< storage::TableHeader> sp_table_header( new storage::TableHeader( function_names));
      storage::Table< double> function_table( sp_table_header);
      function_table.InsertRow( "weights", weights_vector);
      function_table.InsertRow( "value", values_vector);
      function_table.InsertRow( "weighted_value", weighted_values_vector);

      // end
      return function_table;
    }

    //! @brief creates a table with individual function values and their weighted sums for the given argument
    //! the table has scores as rows
    //! @param PROTEIN_MODEL ProteinModel to be used for calculating the function value
    //! @return a table with individual functions and their weighted sums
    storage::Table< double> ProteinModelScoreSum::CreateValueTableVertical
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // cast a pointer to the multiplier data
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
      );

      // set boolean to see if this is multimer
      assemble::ProteinModel multimer_model;
      const bool has_multimers( sp_multiplier.IsDefined());
      if( has_multimers)
      {
        multimer_model = sp_multiplier->operator ()( PROTEIN_MODEL);
      }

      // get the value table
      storage::Table< double> table
      (
        has_multimers ?
          math::SumFunctionMixin< ProteinModel>::CreateValueTableVertical( multimer_model) :
          math::SumFunctionMixin< ProteinModel>::CreateValueTableVertical( PROTEIN_MODEL)
      );

      // end
      return table;
    }

    //! @brief creates a table with individual function values and their weighted sums for the given argument
    //! the table has scores as rows, using readable names instead of scheme
    //! @param PROTEIN_MODEL ProteinModel to be used for calculating the function value
    //! @return a table with individual functions and their weighted sums
    storage::Table< double> ProteinModelScoreSum::CreateSortedReadableTable( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize table header and a table
      util::ShPtr< storage::TableHeader> sp_table_header
      (
        new storage::TableHeader( math::SumFunctionMixin< ProteinModel>::GetValueTableVerticalColumnNames())
      );
      storage::Table< double> function_table( sp_table_header);

      // initialize the sum
      double sum( 0.0);

      // get the combined map
      storage::Map< ProteinModel::TypeEnum, storage::Map< std::string, double> > score_map( CombineSimilarScores());

      // iterate to collect the scores
      for
      (
        const_iterator function_itr( m_Functions.Begin()), function_itr_end( m_Functions.End());
        function_itr != function_itr_end;
        ++function_itr
      )
      {
        // get score
        const util::SiPtr< const ProteinModel> si_score( &*function_itr->Second());
        const ProteinModel &score( *si_score);

        const double this_weighted_value( score( PROTEIN_MODEL) * function_itr->First());

        score_map[ score.GetType()][ score.GetReadableScheme()] += this_weighted_value;

        // update the sums
        sum += this_weighted_value;
      }

      // get the structure based scores
      const storage::Map< ProteinModel::TypeEnum, storage::Map< std::string, double> >::const_iterator struct_find
      (
        score_map.Find( ProteinModel::e_Structure)
      );
      if( struct_find != score_map.End())
      {
        function_table.InsertRow
        (
          ProteinModel::GetTypeName( ProteinModel::e_Structure),
          storage::Vector< double>::Create( util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>())
        );
        for
        (
          storage::Map< std::string, double>::const_iterator score_itr( struct_find->second.Begin()),
            score_itr_end( struct_find->second.End());
          score_itr != score_itr_end; ++score_itr
        )
        {
          function_table.InsertRow
          (
            storage::Table< double>::s_Indentation + score_itr->first,
            storage::Vector< double>::Create( util::GetUndefined< double>(), util::GetUndefined< double>(), score_itr->second)
          );
        }
      }

      // get the sequence based scores
      const storage::Map< ProteinModel::TypeEnum, storage::Map< std::string, double> >::const_iterator seq_find
      (
        score_map.Find( ProteinModel::e_Sequence)
      );
      if( struct_find != score_map.End())
      {
        function_table.InsertRow
        (
          ProteinModel::GetTypeName( ProteinModel::e_Sequence),
          storage::Vector< double>::Create( util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>())
        );
        for
        (
          storage::Map< std::string, double>::const_iterator score_itr( seq_find->second.Begin()),
            score_itr_end( seq_find->second.End());
          score_itr != score_itr_end; ++score_itr
        )
        {
          function_table.InsertRow
          (
            storage::Table< double>::s_Indentation + score_itr->first,
            storage::Vector< double>::Create( util::GetUndefined< double>(), util::GetUndefined< double>(), score_itr->second)
          );
        }
      }

      // add the sum row
      function_table.InsertRow( "Score (total)", storage::Vector< double>::Create( sum, sum, sum));

      // end
      return function_table;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief initialize the score sum using the given score weight map
    //! @param SCORE_WEIGHT_MAP map of scores and corresponding weights
    void ProteinModelScoreSum::AddScoresWithWeights
    (
      const storage::Map< util::ShPtr< ProteinModel>, double> &SCORE_WEIGHT_MAP
    )
    {
      // iterate over weights map provided
      for
      (
        storage::Map< util::ShPtr< ProteinModel>, double>::const_iterator
          map_itr( SCORE_WEIGHT_MAP.Begin()), map_itr_end( SCORE_WEIGHT_MAP.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // if the weight is equal to 0
        if( map_itr->second == double( 0.0))
        {
          // skip this one
          BCL_MessageStd
          (
            "Weight is equal to 0, therefore not adding the following score " + map_itr->first->GetScheme()
          );
        }
        // otherwise
        else
        {
          // add the score
          NewOperand( *map_itr->first, map_itr->second);
        }
      }
    }

    //! @brief initialize the score sum using the given score map and weight map
    //! @param SCORE_MAP map of scoring functions to be used
    //! @param WEIGHT_SET map of function schemes and corresponding weights
    void ProteinModelScoreSum::AddScoresWithWeights
    (
      const storage::Map< std::string, util::ShPtr< ProteinModel> > &SCORE_MAP,
      const storage::Map< std::string, double> &WEIGHT_SET
    )
    {
      // iterate over weights provided
      for
      (
        storage::Map< std::string, double>::const_iterator weight_itr( WEIGHT_SET.Begin()),
          weight_itr_end( WEIGHT_SET.End());
        weight_itr != weight_itr_end;
        ++weight_itr
      )
      {
        // if the weight is equal to 0
        if( weight_itr->second == double( 0.0))
        {
          // skip this one
          BCL_MessageStd
          (
            "Weight is equal to 0, therefore not adding the following score " + weight_itr->first
          );
          continue;
        }

        // look in the score map for this function
        storage::Map< std::string, util::ShPtr< ProteinModel> >::const_iterator
          sp_score_itr( SCORE_MAP.Find( weight_itr->first));

        // make sure it was found
        BCL_Assert( sp_score_itr != SCORE_MAP.End(), "Score was not found in the map " + weight_itr->first);

        // add the score
        NewOperand( *sp_score_itr->second, weight_itr->second);
      }
    }

    //! @brief calculate the score for the given ProteinModel
    //! @brief PROTEIN_MODEL ProteinModel to be evaluated
    //! @return the score for the given model
    double ProteinModelScoreSum::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // cast a pointer to the multiplier data
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
      );

      // set boolean to see if this is multimer
      const bool has_multimers( sp_multiplier.IsDefined());

      // initialize double sum
      double sum( 0.0);

      // calculate the sum function result
      if( has_multimers)
      {
        sum =
         math::SumFunctionMixin< ProteinModel>::operator()( sp_multiplier->operator ()( PROTEIN_MODEL));
      }
      // otherwise score just the protein model
      else
      {
        sum = math::SumFunctionMixin< ProteinModel>::operator()( PROTEIN_MODEL);
      }

      // return sum f(x) = y
      return sum;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write detailed scheme and values to OSTREAM
    //! @param PROTEIN_MODEL ProteinModel to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &ProteinModelScoreSum::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      std::ostream &OSTREAM
    ) const
    {

      // cast a pointer to the multiplier data
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
      );

      // set boolean to see if this is multimer
      const bool has_multimers( sp_multiplier.IsDefined());
      const size_t nr_multimers( has_multimers ? sp_multiplier->GetNumberMultimers() : 1);

      OSTREAM << "has_multimers " << has_multimers << '\n';
      OSTREAM << "nr_multimers " << nr_multimers << '\n';

      for
      (
        const_iterator pair_itr( m_Functions.Begin()), pair_itr_end( m_Functions.End());
        pair_itr != pair_itr_end;
        ++pair_itr
      )
      {
        OSTREAM << pair_itr->Second()->GetScheme() << '\t'
                << nr_multimers * pair_itr->First() * pair_itr->Second()->operator()( PROTEIN_MODEL) << std::endl;

        pair_itr->Second()->WriteDetailedSchemeAndValues( PROTEIN_MODEL, OSTREAM) << '\n';
      }

      //end
      return OSTREAM;
    }

    //! @brief combines scores with the same readable scheme
    //! @return scores with the same readable scheme
    storage::Map< ProteinModel::TypeEnum, storage::Map< std::string, double> > ProteinModelScoreSum::CombineSimilarScores() const
    {
      // initialize map
      storage::Map< ProteinModel::TypeEnum, storage::Map< std::string, double> > score_map;

      // iterate to collect the scores
      for
      (
        const_iterator function_itr( m_Functions.Begin()), function_itr_end( m_Functions.End());
        function_itr != function_itr_end;
        ++function_itr
      )
      {
        // get score
        const util::SiPtr< const ProteinModel> si_score( &*function_itr->Second());
        const ProteinModel &score( *si_score);

        // insert into map
        score_map[ score.GetType()][ score.GetReadableScheme()] = 0;
      }

      // end
      return score_map;
    }

  } // namespace score
} // namespace bcl

