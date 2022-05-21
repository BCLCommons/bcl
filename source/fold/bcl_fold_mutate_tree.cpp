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
#include "fold/bcl_fold_mutate_tree.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_mutates.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    //! @brief conversion to a string from a MutateTypes
    //! @param MUTATE_TYPE the mutate type to get a string for
    //! @return a string representing that mutate type
    const std::string &MutateTree::GetMutateTypeName( const MutateTypes &MUTATE_TYPE)
    {
       static const std::string s_descriptors[] =
       {
         "add", "remove", "swap", "sse", "helix", "strand", "ssepair", "helixpair",
         "helixdomain", "sheet", "domain", "model", "ensemble", "undefined"
       };
       return s_descriptors[ size_t( MUTATE_TYPE)];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateTree::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateTree())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateTree::MutateTree() :
      m_MutateTypeProbabilities(),
      m_MutateProbabilities()
    {
    }

    //! @brief constructor from a weight table
    //! @param WEIGHT_TABLE table with mutate weights
    MutateTree::MutateTree( const storage::Table< double> &WEIGHT_TABLE) :
      m_MutateTypeProbabilities(),
      m_MutateProbabilities()
    {
      // initialize from the table
      InitializeFromTable( WEIGHT_TABLE);
    }

    //! @brief Clone function
    //! @return pointer to new MutateTree
    MutateTree *MutateTree::Clone() const
    {
      return new MutateTree( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateTree::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateTree::GetAlias() const
    {
      static const std::string s_name( "MutateTree");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset all the information in the tree
    void MutateTree::Reset()
    {
      // iterate over the mutate type probabilities and set to 0
      for
      (
        storage::Map< MutateTypesEnum, double>::iterator
          itr( m_MutateTypeProbabilities.Begin()), itr_end( m_MutateTypeProbabilities.End());
        itr != itr_end; ++itr
      )
      {
        itr->second = 0.0;
      }

      // iterate over the mutate probabilities and set to 0
      for
      (
        storage::Map< MutateTypesEnum, storage::Map< Mutate, double> >::iterator
          itr_a( m_MutateProbabilities.Begin()), itr_a_end( m_MutateProbabilities.End());
        itr_a != itr_a_end; ++itr_a
      )
      {
        // iterate over all the mutates for this type
        for
        (
          storage::Map< Mutate, double>::iterator itr_b( itr_a->second.Begin()), itr_b_end( itr_a->second.End());
          itr_b != itr_b_end; ++itr_b
        )
        {
          itr_b->second = 0.0;
        }
      }
    }

    //! @brief sets the probability for a given mutate enum
    //! @param MUTATE_ENUM enum of the mutate of interest
    //! @param PROBABILITY probability of the mutate to be assigned
    void MutateTree::SetMutateProbability( const Mutate &MUTATE_ENUM, const double PROBABILITY)
    {
      // first get the corresponding mutate type
      const MutateTypes mutate_type( ParseMutateType( MUTATE_ENUM.GetName()));

      // make sure it is not undefined
      BCL_Assert
      (
        mutate_type != s_NumberMutateTypes, "The given mutate name does not conform to standards " +
        MUTATE_ENUM.GetName()
      );

      // set the probability
      m_MutateProbabilities[ mutate_type][ MUTATE_ENUM] = PROBABILITY;
    }

    //! @brief sets the default probability for a given mutate enum
    //! @param MUTATE_ENUM enum of the mutate of interest
    //! @param PROBABILITY probability of the mutate to be assigned
    void MutateTree::SetDefaultMutateProbability( const Mutate &MUTATE_ENUM, const double PROBABILITY)
    {
      // first get the corresponding mutate type
      const MutateTypes mutate_type( ParseMutateType( MUTATE_ENUM.GetName()));

      // make sure it is not undefined
      BCL_Assert
      (
        mutate_type != s_NumberMutateTypes, "The given mutate name does not conform to standards " +
        MUTATE_ENUM.GetName()
      );

      // set the probability
      m_MutateProbabilities[ mutate_type][ MUTATE_ENUM] = PROBABILITY;

      // set the set status to not set
      SetMutateSet( MUTATE_ENUM, false);
    }

    //! @brief gets the probability for a given mutate enum
    //! @param MUTATE_ENUM enum of the mutate of interest
    //! @return probability of the mutate
    double MutateTree::GetMutateProbability( const Mutate &MUTATE_ENUM) const
    {
      // first get the corresponding mutate type
      const MutateTypes mutate_type( ParseMutateType( MUTATE_ENUM.GetName()));

      // get an iterator to the mutate type
      storage::Map< MutateTypesEnum, storage::Map< Mutate, double> >::const_iterator type_itr
      (
        m_MutateProbabilities.Find( mutate_type)
      );
      if( type_itr == m_MutateProbabilities.End())
      {
        return 0.0;
      }

      // get an iterator to the mutate
      storage::Map< Mutate, double>::const_iterator mutate_itr( type_itr->second.Find( MUTATE_ENUM));

      // end
      return mutate_itr == type_itr->second.End() ? 0.0 : mutate_itr->second;
    }

    //! @brief sets the probability for a given mutate
    //! @param MUTATE_NAME Name of the mutate of interest
    //! @param PROBABILITY probability of the mutate to be assigned
    void MutateTree::SetMutateProbability( const std::string &MUTATE_NAME, const double PROBABILITY)
    {
      // get the corresponding enum
      const Mutate mutate( GetMutates().GetEnumFromName( MUTATE_NAME));

      BCL_Assert( mutate.IsDefined(), "no enum with name: " + MUTATE_NAME);

      // call with enum
      SetMutateProbability( mutate, PROBABILITY);
    }
    //! @brief get the probability of a given mutate type
    //! @param MUTATE_TYPE mutate type of interest
    //! return probability of mutate type
    double MutateTree::GetMutateTypeProbability( const MutateTypes &MUTATE_TYPE) const
    {
      // get an iterator to the mutate type
      storage::Map< MutateTypesEnum, double>::const_iterator type_itr
        (
         m_MutateTypeProbabilities.Find( MUTATE_TYPE)
         );

      // make sure it was found
      BCL_Assert
        (
         type_itr != m_MutateTypeProbabilities.End(),
         "The given mutate name does not conform to standards "
         );

      // end
      return type_itr->second;
    }

    //! @brief set the probabilities for a specific mutate type
    //! @param MUTATE_TYPE MutateTypes of interest
    //! @param PROBABILITY Probability to set to
    void MutateTree::SetMutateTypeProbability( const MutateTypes &MUTATE_TYPE, const double PROBABILITY)
    {
      m_MutateTypeProbabilities[ MUTATE_TYPE] = PROBABILITY;
      m_MutateTypeSet[ MUTATE_TYPE] = true;
    }

    //! @brief set the default probabilities for a specific mutate type
    //! @param MUTATE_TYPE MutateTypes of interest
    //! @param PROBABILITY Probability to set to
    void MutateTree::SetDefaultMutateTypeProbability( const MutateTypes &MUTATE_TYPE, const double PROBABILITY)
    {
      m_MutateTypeProbabilities[ MUTATE_TYPE] = PROBABILITY;
      m_MutateTypeSet[ MUTATE_TYPE] = false;
    }

    //! @brief gets the set status for a given mutate type
    //! @param MUTATE_TYPE type of the mutate of interest
    //! @return set status of the mutate type
    bool MutateTree::GetMutateTypeSet( const MutateTypes &MUTATE_TYPE) const
    {
      // get an iterator to the mutate type
      storage::Map< MutateTypesEnum, bool>::const_iterator type_itr
        (
         m_MutateTypeSet.Find( MUTATE_TYPE)
         );
      if( type_itr == m_MutateTypeSet.End())
        {
          return false;
        }
      // end
      return type_itr->second;
    }

    //! @brief set the mutate type to having been set
    //! @param MUTATE_TYPE mutate type of interest
    //! @param SET whether or not it has been set
    void MutateTree::SetMutateTypeSet( const MutateTypes &MUTATE_TYPE, const bool SET)
    {
      m_MutateTypeSet[ MUTATE_TYPE] = SET;
    }

    //! @brief get whether or not the probability of the given mutate type is zero
    //! @param MUTATE_TYPE mutate type of interest
    //! @return whether or not the probability is zero
    bool MutateTree::GetMutateTypeZero( const MutateTypes &MUTATE_TYPE) const
    {
      storage::Map< MutateTypesEnum, double>::const_iterator type_itr
      (
        m_MutateTypeProbabilities.Find( MUTATE_TYPE)
      );

      return type_itr == m_MutateTypeProbabilities.End() || type_itr->second == 0.0;
    }

    //! @brief get whether or not a given mutate has been set
    //! @param MUTATE mutate of interest
    //! @return whether or not this mutate has been set
    bool MutateTree::GetMutateSet( const Mutate &MUTATE) const
    {
      // get an iterator to the mutate type
      storage::Map< Mutate, bool>::const_iterator set_itr
      (
        m_MutateSet.Find( MUTATE)
      );

      return set_itr != m_MutateSet.End() && set_itr->second;
    }

    //! @brief set the mutate to having been set
    //! @param MUTATE mutate of interest
    //! @param SET whether or not it has been set
    void MutateTree::SetMutateSet( const Mutate &MUTATE, const bool SET)
    {
      m_MutateSet[ MUTATE] = SET;
    }

    //! @brief get whether or not the probability of the given mutate is zero
    //! @param MUTATE mutate of interest
    //! @return whether or not the probability is zero
    bool MutateTree::GetMutateZero( const Mutate &MUTATE_ENUM) const
    {
      // first get the corresponding mutate type
      const MutateTypes mutate_type( ParseMutateType( MUTATE_ENUM.GetName()));

      // get an iterator to the mutate type
      storage::Map< MutateTypesEnum, storage::Map< Mutate, double> >::const_iterator type_itr
        (
         m_MutateProbabilities.Find( mutate_type)
         );

      // make sure it was found
      if( type_itr == m_MutateProbabilities.End())
      {
        return 0.0;
      }

      // get an iterator to the mutate
      storage::Map< Mutate, double>::const_iterator mutate_itr( type_itr->second.Find( MUTATE_ENUM));

      return mutate_itr == type_itr->second.End() || mutate_itr->second == 0.0;
    }

    //! @brief parses the given mutate name and returns the mutate type
    //! @param MUTATE_NAME name of the mutate
    //! @return corresponding mutate type
    MutateTree::MutateTypes MutateTree::ParseMutateType( const std::string &MUTATE_NAME)
    {
      // get the substring until the first "_"
      const size_t first_underscore( MUTATE_NAME.find_first_of( '_'));

      // if it was not found
      BCL_Assert( first_underscore != std::string::npos, "Following mutate does not exist: " + MUTATE_NAME);

      // form the substring
      const std::string mutate_type_name( MUTATE_NAME.substr( 0, first_underscore));

      // find the matching enum (returns s_NumberMutateTypes if none have this name)
      return MutateTypesEnum( mutate_type_name);
    }

    //! @brief create a table from mutate tree
    //! @return table created from the mutate tree
    storage::Table< double> MutateTree::CreateTable() const
    {
      // initialize table
      storage::Table< double> table( storage::Vector< std::string>::Create( "weights", "final_weights"));

      // iterate over the mutate types
      for
      (
        storage::Map< MutateTypesEnum, double>::const_iterator
          map_itr( m_MutateTypeProbabilities.Begin()), map_itr_end( m_MutateTypeProbabilities.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // insert the name
        table.InsertRow
        (
          GetMutateTypeName( map_itr->first),
          storage::Vector< double>::Create( map_itr->second, map_itr->second)
        );
      }

      // iterate over the mutate tree
      for
      (
        storage::Map< MutateTypesEnum, storage::Map< Mutate, double> >::const_iterator
          node_itr( m_MutateProbabilities.Begin()), node_itr_end( m_MutateProbabilities.End());
        node_itr != node_itr_end; ++node_itr
      )
      {
        // create a reference on this leaf map
        const storage::Map< Mutate, double> &this_map( node_itr->second);

        // now iterate over the leaves
        for
        (
          storage::Map< Mutate, double>::const_iterator
            mutate_itr( this_map.Begin()), mutate_itr_end( this_map.End());
          mutate_itr != mutate_itr_end; ++mutate_itr
        )
        {
          // get the probability of the mutate type
          const storage::Map< MutateTypesEnum, double>::const_iterator type_itr
          (
            m_MutateTypeProbabilities.Find( node_itr->first)
          );
          BCL_Assert
          (
            type_itr != m_MutateTypeProbabilities.End(),
            "Unable to find mutate type for " + mutate_itr->first.GetName()
          );

          // insert the leaf into the table
          table.InsertRow
          (
            mutate_itr->first.GetName(),
            storage::Vector< double>::Create
            (
              mutate_itr->second,
              type_itr->second * mutate_itr->second
            )
          );
        }
      }

      // end
      return table;
    }

    //! @brief creates the mutate function using all the mutates in the mutate tree
    //! @return the mutate function constructed from all the mutates in the mutate tree
    util::ShPtr< math::MutateInterface< assemble::ProteinModel> > MutateTree::ConstructMutates() const
    {
      // initialize object distribution for all nodes
      math::ObjectProbabilityDistribution< math::MutateInterface< assemble::ProteinModel> > all_nodes;

      // iterate over all nodes
      for
      (
        storage::Map< MutateTypesEnum, double>::const_iterator
          map_itr( m_MutateTypeProbabilities.Begin()), map_itr_end( m_MutateTypeProbabilities.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // create variables to hold the type and the probability
        const MutateTypes mutate_type( map_itr->first);
        const std::string mutate_type_name( GetMutateTypeName( mutate_type));
        const double mutate_type_prob( map_itr->second);

        // see if there is an entry for this mutate type
        storage::Map< MutateTypesEnum, storage::Map< Mutate, double> >::const_iterator itr( m_MutateProbabilities.Find( mutate_type));

        // if not found
        if( itr == m_MutateProbabilities.End())
        {
          // then skip to next one
          continue;
        }

        // create reference on the mutate map for this type
        const storage::Map< Mutate, double> &this_mutate_map( itr->second);

        BCL_MessageVrb
        (
          "mutate_type: " + mutate_type_name + " " + util::Format()( mutate_type_prob)
        );

        // warn user if node probability is set to 0 but there are corresponding leaves
        if( mutate_type_prob == 0.0)
        {
          if( !this_mutate_map.IsEmpty())
          {
            BCL_MessageVrb
            (
              "The mutate node " + mutate_type_name +
              " has a 0 probability, therefore following mutates will not be used\n" + util::Format()( this_mutate_map)
            );
          }
        }
        else if( this_mutate_map.IsEmpty())
        {
          BCL_MessageVrb( "No mutates in this node therefore skipping " + mutate_type_name);
        }
        else
        {
          // build a object probability distribution for this node
          math::ObjectProbabilityDistribution< math::MutateInterface< assemble::ProteinModel> > this_node;

          // iterate over mutates of this type and add them to this node
          for
          (
            storage::Map< Mutate, double>::const_iterator
              mutate_itr( this_mutate_map.Begin()), mutate_itr_end( this_mutate_map.End());
            mutate_itr != mutate_itr_end; ++mutate_itr
          )
          {
            // store the mutate name
            const std::string mutate_name( mutate_itr->first.GetName());
            const double mutate_prob( mutate_itr->second);

            BCL_MessageVrb
            (
              "mutate " + mutate_name + " " + util::Format()( mutate_prob)
            );

            // if 0.0 probability skip over
            if( mutate_prob == 0.0)
            {
              BCL_MessageVrb( "Skipping this mutate with 0 probability " + mutate_name);
              continue;
            }

            // insert into the map the mutate and the probability
            this_node.PushBack( mutate_prob, **mutate_itr->first);
          }
          // construct decision node push it into all_nodes
          all_nodes.PushBack( mutate_type_prob, math::MutateDecisionNode< assemble::ProteinModel>( this_node, map_itr->first.GetString()));
        }
      }

      // construct the decision node from all nodes and return
      return util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
      (
        new math::MutateDecisionNode< assemble::ProteinModel>( all_nodes, "MutateDecisionNode")
      );
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateTree::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Stores mutates with associated probabilities for MCM approximation");
      serializer.AddInitializer
      (
        "mutate type probabilities",
        "assignment of probabilities to the different mutate types",
        io::Serialization::GetAgent( &m_MutateTypeProbabilities)
      );
      serializer.AddInitializer
      (
        "mutate probabilities",
        "assignment of probabilities to the different mutates",
        io::Serialization::GetAgent( &m_MutateProbabilities)
      );

      return serializer;
    }

    //! @brief merge two mutate trees
    void MutateTree::Merge( const MutateTree &MUTATE_TREE)
    {
      // set mutate type probabilities
      storage::Map< MutateTypesEnum, double> mutate_types( MUTATE_TREE.GetMutateTypeProbabilities());
      for( auto type_itr( mutate_types.Begin()), type_itr_end( mutate_types.End()); type_itr != type_itr_end; ++type_itr)
        {
          // if it has not been set or the new tree is setting it to zero
          if( !GetMutateTypeSet( type_itr->first) || MUTATE_TREE.GetMutateTypeZero( type_itr->first))
            {
              // set the tree's probability to that of MUTATE_TREE
              const double type_prob( MUTATE_TREE.GetMutateTypeProbability( type_itr->first));
              SetMutateTypeProbability( type_itr->first, type_prob);
              SetMutateTypeSet( type_itr->first, true);
            }
        }

      // set mutate probabilities according to the second tree
      storage::Map< MutateTypesEnum, storage::Map< Mutate, double> > mutates( MUTATE_TREE.GetMutateProbabilities());
      for( auto node_itr( mutates.Begin()), node_itr_end( mutates.End()); node_itr != node_itr_end; node_itr++)
      {
        // create a reference on this leaf map
         const storage::Map< Mutate, double> &this_map( node_itr->second);

        //iterate over the leaves
         for( auto mutate_itr( this_map.Begin()), mutate_itr_end( this_map.End()); mutate_itr != mutate_itr_end; ++mutate_itr)
         {
           // set new mutate probabilities if they haven't already been set or the new probability is zero
           if( !GetMutateSet( mutate_itr->first) || MUTATE_TREE.GetMutateZero( mutate_itr->first))
           {
             const double mutate_prob( MUTATE_TREE.GetMutateProbability( mutate_itr->first));
             SetMutateProbability( mutate_itr->first, mutate_prob);
             SetMutateSet( mutate_itr->first, true);
           }
         }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initializes the mutate tree using a given table
    //! @param WEIGHT_TABLE table that contains the weights for mutates
    void MutateTree::InitializeFromTable( const storage::Table< double> &WEIGHT_TABLE)
    {
      // if the weight table does not have a row named weights
      if( !WEIGHT_TABLE.HasRow( "weights"))
      {
        // make sure at least it has column weights
        BCL_Assert
        (
          WEIGHT_TABLE.GetHeader().HasColumn( "weights"),
          "The given mutate weight table has no rows or columns named \"weights\""
        );

        // then transpose the table
        storage::Table< double> weight_table_transposed( WEIGHT_TABLE.GetTransposedTable());

        // now call this function with the transposed table
        return InitializeFromTable( weight_table_transposed);
      }

      // get the map from the table
      storage::Map< std::string, double> weights_read( WEIGHT_TABLE[ "weights"].ConvertToMap());

      // iterate over each mutate type
      for( size_t i( 0); i < s_NumberMutateTypes; ++i)
      {
        // construct mutate type
        MutateTypes this_type = MutateTypes( i);
        MutateTypesEnum wrapper( this_type);

        // try to find the probability for this mutate type
        storage::Map< std::string, double>::iterator type_itr( weights_read.Find( wrapper));

        // if there was no probability
        if( type_itr == weights_read.End())
        {
          // insert with 0 probability
          m_MutateTypeProbabilities[ this_type] = 0.0;
        }
        // otherwise
        else
        {
          // insert the weight
          m_MutateTypeProbabilities[ this_type] = type_itr->second;

          // remove from the list
          weights_read.RemoveElement( type_itr);
        }
      }

      // now iterate over remaining ones - they should be all mutates
      for
      (
        storage::Map< std::string, double>::const_iterator
          map_itr( weights_read.Begin()), map_itr_end( weights_read.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // store the string and probability
        const std::string mutate_name( map_itr->first);
        const double mutate_prob( map_itr->second);

        // parse the mutate type and type name
        const MutateTypes this_type( ParseMutateType( mutate_name));
        const std::string this_typename( GetMutateTypeName( this_type));

        // make sure it is not an undefined mutate type
        BCL_Assert( this_type != s_NumberMutateTypes, "this mutate name is not valid! " + mutate_name);

        // find the corresponding mutate
        const Mutate &mutate( GetMutates().GetEnumFromName( mutate_name));

        // make sure it is not undefined
        BCL_Assert( mutate != GetMutates().e_Undefined, "There is no mutate with the given name " + mutate_name);

        // if the node has zero probability
        if( m_MutateTypeProbabilities[ this_type] == 0.0)
        {
          // warn user
          BCL_MessageVrb
          (
            "Following mutate won't be used since the node it belongs has a 0 probability " + mutate_name
          );
        }
        // if the node has non zero probability then insert
        else
        {
          // otherwise insert into corresponding node
          m_MutateProbabilities[ this_type][ mutate] = mutate_prob;
        }
      }
    }

  } // namespace fold

} // namespace bcl
