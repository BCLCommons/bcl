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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "fold/bcl_fold_default_flags.h"
#include "io/bcl_io_directory_entry.h"
#include "math/bcl_math_mutate_combine.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "math/bcl_math_sum_function.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_handler_atom_distance_assigned.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_add.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_aa_type.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_coordinate_exclusion.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_euclidian_distance.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_exposure.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_sse_size.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_triangulation.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_remove.h"
#include "score/bcl_score_data_set_pairwise_bipolar.h"
#include "score/bcl_score_data_set_pairwise_coordinate_exclusion.h"
#include "score/bcl_score_data_set_pairwise_coordinate_triangulation.h"
#include "score/bcl_score_data_set_pairwise_data_density.h"
#include "score/bcl_score_data_set_pairwise_distance_change_magnitude.h"
#include "score/bcl_score_data_set_pairwise_euclidian_distance.h"
#include "score/bcl_score_data_set_pairwise_residue_type_exclusion.h"
#include "score/bcl_score_data_set_pairwise_sequence_separation.h"
#include "score/bcl_score_data_set_pairwise_size.h"
#include "score/bcl_score_data_set_pairwise_sse_center.h"
#include "score/bcl_score_data_set_pairwise_sse_connection.h"
#include "score/bcl_score_data_set_pairwise_sse_size.h"
#include "score/bcl_score_data_set_pairwise_sse_term.h"
#include "score/bcl_score_data_set_pairwise_structural_exposure.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{

  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OptimizeDataSetPairwise
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @author alexanns
    //! @date May 12, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API OptimizeDataSetPairwise :
      public Interface
    {

    private:

      typedef util::ShPtr< math::MutateInterface< restraint::DataSetPairwise> > Mutate;

      //! input file with distances between objects
      util::ShPtr< command::FlagInterface> m_WriteRestraintsBCLFormat;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class MathFunctions
      //! @brief TODO: add a brief comment
      //! @details TODO: add an detailed description to this class
      //!
      //! @author alexanns
      //! @date May 12, 2011
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      template< typename t_Argument, typename t_Return>
      class MathFunctions :
        public util::Enumerate< util::ShPtr< math::FunctionInterfaceSerializable< t_Argument, t_Return> >, MathFunctions< t_Argument, t_Return> >
      {
        friend class util::Enumerate< util::ShPtr< math::FunctionInterfaceSerializable< t_Argument, t_Return> >, MathFunctions< t_Argument, t_Return> >;

      private:

        //! @brief default constructor
        MathFunctions()
        {
        }

      //////////
      // data //
      //////////

      public:

        typedef typename util::Enumerate< util::ShPtr< math::FunctionInterfaceSerializable< t_Argument, t_Return> >, MathFunctions< t_Argument, t_Return> >::const_iterator const_iterator;

        //! single instance of that class
        static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name
        //! @return the class name as const ref std::string
        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName( *this);
        }

      ////////////////
      // operations //
      ////////////////

        //! @brief adds a score to the enumerated scores
        //! @param SCORE Score to add
        typename MathFunctions< t_Argument, t_Return>::EnumType AddFunction( const util::ShPtr< math::FunctionInterfaceSerializable< t_Argument, t_Return> > &FUNCTION)
        {
          return AddEnum( FUNCTION->GetScheme(), FUNCTION);
        }

      private:

        //! @brief function for adding a new enum
        //! @param NAME name of the current enum
        //! @param OBJECT object to be enumerated
        typename MathFunctions< t_Argument, t_Return>::EnumType &AddEnum
        (
          const std::string &NAME,
          const util::ShPtr< math::FunctionInterfaceSerializable< t_Argument, t_Return> > &FUNCTION
        )
        {
          bool name_exists( util::Enumerate< util::ShPtr< math::FunctionInterfaceSerializable< t_Argument, t_Return> >, MathFunctions< t_Argument, t_Return> >::HaveEnumWithName( NAME));
          // make sure an enum with the given name does not exist already
          BCL_Assert( !name_exists, "A score enum with the given name already exists: " + NAME);

          // call the add enum
          return util::Enumerate< util::ShPtr< math::FunctionInterfaceSerializable< t_Argument, t_Return> >, MathFunctions< t_Argument, t_Return> >::AddEnum( NAME, FUNCTION);
        }

      ///////////////
      // operators //
      ///////////////

      //////////////////////
      // input and output //
      //////////////////////

      //////////////////////
      // helper functions //
      //////////////////////

      }; // class Scores

      //! @brief construct on access function for all Scores
      //! @return reference to only instances of Scores
      template< typename t_Argument, typename t_Return> static
      MathFunctions< t_Argument, t_Return> &GetMathFunctions()
      {
        return MathFunctions< t_Argument, t_Return>::GetEnums();
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class MathFunctionsWeightSet
      //! @brief TODO: add a brief comment
      //! @details TODO: add an detailed description to this class
      //!
      //! @author alexanns
      //! @date May 15, 2011
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template< typename t_Argument, typename t_Return>
      class MathFunctionsWeightSet :
        public util::ObjectInterface
      {

      private:

        //! @typedef for enum
        typedef util::Enum< util::ShPtr< math::FunctionInterfaceSerializable< t_Argument, t_Return> >, MathFunctions< t_Argument, t_Return> > MathFunction;

      //////////
      // data //
      //////////

        //! list of enum and associated weights
        storage::List< storage::Pair< MathFunction, double> > m_WeightMap;

      public:

        //! single instance of that class
        static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief default constructor
        MathFunctionsWeightSet() :
          m_WeightMap()
        {
        }

        //! @brief constructor from a table
        //! @param TABLE that contains the weights
        MathFunctionsWeightSet( const MathFunctions< t_Argument, t_Return> &MATH_FUNCTIONS, const double DEFAULT_WEIGHT) :
          m_WeightMap()
        {
          InitializeFromMathFunctions( MATH_FUNCTIONS, DEFAULT_WEIGHT);
        }

        //! @brief constructor from a table
        //! @param TABLE that contains the weights
        MathFunctionsWeightSet( const storage::Table< double> &TABLE) :
          m_WeightMap()
        {
          InitializeFromTable( TABLE);
        }

        //! @brief Clone function
        //! @return pointer to new WeightSet
        MathFunctionsWeightSet *Clone() const
        {
          return new MathFunctionsWeightSet( *this);
        }

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name
        //! @return the class name as const ref std::string
        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName( *this);
        }

      ////////////////
      // operations //
      ////////////////

        //! @brief get weight for a given score
        //! @param FUNCTION Score of interest
        //! @return weight for the given FUNCTION
        double GetWeight( const MathFunction &FUNCTION)
        {
          typename storage::List< storage::Pair< MathFunction, double> >::const_iterator itr( std::find( m_WeightMap.Begin(), m_WeightMap.End(), FUNCTION));

          return itr == m_WeightMap.End() ? util::GetUndefinedDouble() : itr->Second();
        }

        //! @brief sets the weight of the given score
        //! @param FUNCTION_NAME name of the function
        //! @param WEIGHT Weight to be asssigned to function
        void SetWeight
        (
          const std::string &FUNCTION_NAME,
          const double WEIGHT
        )
        {
          // find the enum
          const MathFunction &function( GetMathFunctions< t_Argument, t_Return>().GetEnumFromName( FUNCTION_NAME));

          // make sure the a score with the given name exists
          MathFunction undefined( GetMathFunctions< t_Argument, t_Return>().e_Undefined);
          BCL_Assert( function != undefined, "There is no score with given name " + FUNCTION_NAME);

          typename storage::List< storage::Pair< MathFunction, double> >::iterator itr
          (
            std::find( m_WeightMap.Begin(), m_WeightMap.End(), function)
          );

          if( itr != m_WeightMap.End())
          {
            // update the weight
            itr->Second() = WEIGHT;
          }
        }

        //! @brief constructs the scores and returns it
        //! @return constructed score
        util::ShPtr< math::SumFunction< t_Argument, t_Return> > ConstructSumFunction() const
        {
          util::ShPtr< math::SumFunction< t_Argument, t_Return> > sum_function
          (
            new math::SumFunction< t_Argument, t_Return>()
          );

          // iterate over weights map
          for
          (
            typename storage::List< storage::Pair< MathFunction, double> >::const_iterator
              map_itr( m_WeightMap.Begin()), map_itr_end( m_WeightMap.End());
            map_itr != map_itr_end; ++map_itr
          )
          {
            // if the weight is equal to 0
            if( map_itr->Second() == double( 0.0))
            {
              // skip this one
              BCL_MessageStd
              (
                "Weight is equal to 0, therefore not adding the following function " + map_itr->First().GetName()
              );
            }
            // otherwise
            else
            {
              // add the score
              sum_function->NewOperand( **map_itr->First(), map_itr->Second());
            }
          }

          return sum_function;
        }

        math::ObjectProbabilityDistribution< math::MutateInterface< t_Argument> >
        ConstructObjectProbabilityDistribution() const
        {
          typedef math::ObjectProbabilityDistribution< math::MutateInterface< t_Argument> > t_Object;
          typedef typename t_Object::Assignment Assignment;
          storage::Vector< Assignment> functions;

          // iterate over weights map
          for
          (
            typename storage::List< storage::Pair< MathFunction, double> >::const_iterator
              map_itr( m_WeightMap.Begin()), map_itr_end( m_WeightMap.End());
            map_itr != map_itr_end; ++map_itr
          )
          {
            const math::MutateInterface< t_Argument> &function_tmp
            (
              dynamic_cast< const math::MutateInterface< t_Argument> &>( **map_itr->First())
            );
            functions.PushBack( Assignment( map_itr->Second(), function_tmp));
          }

          return functions;
        }

        util::ShPtr< math::MutateCombine< t_Argument> > ConstructMutateCombine() const
        {
          util::ShPtrList< math::FunctionInterfaceSerializable< t_Argument, math::MutateResult< t_Argument> > > mutates;

          // iterate over weights map
          for
          (
            typename storage::List< storage::Pair< MathFunction, double> >::const_iterator
              map_itr( m_WeightMap.Begin()), map_itr_end( m_WeightMap.End());
            map_itr != map_itr_end; ++map_itr
          )
          {
            // make ShPtrList with necessary number of mutates and append it to mutates
            const util::ShPtrList< math::FunctionInterfaceSerializable< t_Argument, math::MutateResult< t_Argument> > >
              current_mutates( map_itr->Second(), map_itr->First());

            mutates.Append( current_mutates);
          }

          util::ShPtr< math::MutateCombine< t_Argument> > mutate_combine
          (
            new math::MutateCombine< t_Argument>( mutates, false, math::MutateCombine< t_Argument>::GetDefaultScheme())
          );

          return mutate_combine;
        }

        //! @brief writes the enum schemes in table format with weights of 0
        //! @return ostream to which the enum is written
        storage::Table< double> CreateTable() const
        {
          // create a vector to hold column names
          storage::Vector< std::string> column_names;
          // create a vector to hold the weightst
          storage::Vector< double> weights;

          // iterate over the map
          for
          (
            typename storage::List< storage::Pair< MathFunction, double> >::const_iterator
              map_itr( m_WeightMap.Begin()), map_itr_end( m_WeightMap.End());
            map_itr != map_itr_end; ++map_itr
          )
          {
            // insert the name and the weight
            column_names.PushBack( map_itr->First().GetName());
            weights.PushBack( map_itr->Second());
          }

          // create a table
          storage::Table< double> table( column_names);
          // insert a row
          table.InsertRow( "weights", weights);

          // end
          return table;
        }

        //! @brief set all weights to given weight
        //! @param WEIGHT weight to which all weights will be set
        void SetWeights( const double WEIGHT)
        {
          // iterate over the map
          for
          (
            typename storage::List< storage::Pair< MathFunction, double> >::iterator map_itr( m_WeightMap.Begin()), map_itr_end( m_WeightMap.End());
            map_itr != map_itr_end; ++map_itr
          )
          {
            // set the weight to WEIGHT
            map_itr->Second() = WEIGHT;
          }
        }

        //! @brief resets all weights to zero
        void Reset()
        {
          SetWeights( 0);
        }

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
        std::istream &Read( std::istream &ISTREAM)
        {
          // read members
          io::Serialize::Read( m_WeightMap, ISTREAM);

          // return the stream
          return ISTREAM;
        }

        //! @brief write to std::ostream
        //! @param OSTREAM outputstream to write to
        //! @param INDENT number of indentations
        //! @return outputstream which was written to
        std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
        {
          // write members
          io::Serialize::Write( m_WeightMap, OSTREAM);

          // return the stream
          return OSTREAM;
        }

      //////////////////////
      // helper functions //
      //////////////////////

      private:

        //! @brief initialize from enumerator
        //! @param MATH_FUNCTIONS enumerator to initialize from
        void InitializeFromMathFunctions( const MathFunctions< t_Argument, t_Return> &MATH_FUNCTIONS, const double DEFAULT_WEIGHT)
        {
          // iterate through the MATH_FUNCTIONS
          for
          (
            typename MathFunctions< t_Argument, t_Return>::const_iterator
              function_itr( MATH_FUNCTIONS.Begin()), function_itr_end( MATH_FUNCTIONS.End());
            function_itr != function_itr_end;
            ++function_itr
          )
          {
            m_WeightMap.PushBack( storage::Pair< MathFunction, double>( *function_itr, DEFAULT_WEIGHT));
          }
        }

        //! @brief initialize this ScoreWeightSet from a table
        //! @param TABLE Table that contains the weightset
        void InitializeFromTable( const storage::Table< double> &TABLE)
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
          storage::List< storage::Pair< std::string, double> > string_map( TABLE[ "weights"].ConvertToPairList());

          // iterate over the map
          for
          (
            typename storage::List< storage::Pair< std::string, double> >::const_iterator
              map_itr( string_map.Begin()), map_itr_end( string_map.End());
            map_itr != map_itr_end; ++map_itr
          )
          {
            BCL_MessageStd( "To weightset adding " + map_itr->First());

            // get the score enum with the corresponding name
            const MathFunction &function( GetMathFunctions< t_Argument, t_Return>().GetEnumFromName( map_itr->First()));

            // make sure there is a score with the given string
            MathFunction undefined( GetMathFunctions< t_Argument, t_Return>().e_Undefined);
            BCL_Assert( function != undefined, "There is no score with the name " + map_itr->First());

            // set the weight for this score in the weight map
            m_WeightMap.PushBack( storage::Pair< MathFunction, double>( function, map_itr->Second()));
          }

          BCL_MessageStd
          (
            "After initializing from table weightset size is "
            + util::Format()( m_WeightMap.GetSize())
          );
        }

      }; // class MathFunctionsWeightSet

    //////////
    // data //
    //////////

      //! @brief return command line flag for specifying the residues that should not be considered
      //! @return command line flag for specifying the residues that should not be considered
      static util::ShPtr< command::FlagInterface> &GetFlagAATypeExclusions();

      //! @brief return command line flag for specifying the desired possible size range of the data set
      //! @return command line flag for specifying the desired possible size range of the data set
      static util::ShPtr< command::FlagStatic> &GetFlagDataSetSizeRange();

      //! @brief return command line flag for specifying files containing coordinates for exclusion
      //! @return command line flag for specifying coordinates for exclusion
      static util::ShPtr< command::FlagDynamic> &GetFlagCoordinateExclusionFiles();

      //! @brief return command line flag for specifying the x,y,z columns containing coordinates for exclusion
      //! @return command line flag for specifying the x,y,z columns with coordinates for exclusion
      static util::ShPtr< command::FlagStatic> &GetFlagCoordinateExclusionFileColumns();

      //! @brief specify the radius around exclusion coordinates that will be used to exclude residues
      //! @return command line flag for specifying the radius around exclusion coordinates that will be used
      static util::ShPtr< command::FlagStatic> &GetFlagCoordinateExclusionRadius();

      //! @brief return command line flag for specifying files containing lists of pdbs to make ensembles
      //! @return command line flag for specifying files containing lists of pdbs to make ensembles
      static util::ShPtr< command::FlagDynamic> &GetFlagEnsembleFiles();

      //! @brief specify the radius around residues that will keep close-by residues from being considered
      //! @return command line flag for the radius around residues to keep close-by residues from being considered
      static util::ShPtr< command::FlagStatic> &GetFlagTriangulationRadius();

      //! @brief specify the min and max distances that are desired to be possible for the dataset
      //! @return command line flag for specifying the min and max distances
      static util::ShPtr< command::FlagStatic> &GetFlagEuclidianDistance();

      //! @brief return command line flag for writing scores to a file
      //! @return command line flag for writing scores to a file
      static util::ShPtr< command::FlagInterface> &GetFlagWritePossibleScores();

      //! @brief return command line flag for specifying desired residue exposure
      //! @return command line flag for specifying desired residue exposure
      static util::ShPtr< command::FlagInterface> &GetFlagExposure();

      //! @brief return command line flag for writing mutates to a file
      //! @return command line flag for writing mutates to a file
      static util::ShPtr< command::FlagInterface> &GetFlagWritePossibleMutates();

      //! @brief return command line flag for reading mutates from a table file for creating starting dataset
      //! @return command line flag for reading mutates from a table file for creating starting dataset
      static util::ShPtr< command::FlagInterface> &GetFlagReadMutatesWeightsStart();

      //! @brief command line flag for reading scores from a table file to use for optimizing the starting data set
      //! @return command line flag for reading mutates from a table file to use for optimizing the starting data set
      static util::ShPtr< command::FlagInterface> &GetFlagReadScoresWeightsOptimization();

      //! @brief return command line flag for reading mutates from a table file for optimizing the starting data set
      //! @return command line flag for reading mutates from a table file for optimizing the starting data set
      static util::ShPtr< command::FlagInterface> &GetFlagReadMutatesWeightsOptimization();

      //! @brief return command line flag for reading mutates from a table file for filtering optimized data set
      //! @return command line flag for reading mutates from a table file for filtering optimized data set
      static util::ShPtr< command::FlagInterface> &GetFlagReadMutatesWeightsEnd();

      //! @brief return command line flag for writing a pymol script to show the distances
      //! @return command line flag for writing a pymol script to show the distances
      static util::ShPtr< command::FlagInterface> &GetFlagPymolOutput();

      //! @brief the filename containing pdbs of the structures that should be used to calculate the distance restraints
      //! @return command line flag for specifying pdbs to be used for calculated restraint distances
      static util::ShPtr< command::FlagInterface> &GetFlagRestraintDistanceEnsemble();

      //! @brief flag for specifying the number of desired restraints as fraction of number of residues in sses
      //! @return flag for specifying the number of desired restraints as fraction of number of residues in sses
      static const util::ShPtr< command::FlagInterface> &GetFlagDataSetSizeFractionOfPool();

      //! @brief flag for specifying the number of desired restraints as connecting all sses a given number of times
      //! @return flag for specifying the number of desired restraints as connecting all sses a given number of times
      static const util::ShPtr< command::FlagInterface> &GetFlagDataSetSizeConnectSSEsInPool();

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      OptimizeDataSetPairwise();

      //! @brief Clone function
      //! @return pointer to new OptimizeDataSetPairwise
      OptimizeDataSetPairwise *Clone() const
      {
        return new OptimizeDataSetPairwise( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

    ////////////////
    // operations //
    ////////////////

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
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief initializes and fills the enumerator with all scores
      void InitializeScores() const;

      //! @brief initializes and fills the enumerator with all mutates
      void InitializeMutates() const;

      //! @brief initialize mutates that will not change for a second ensemble
      void InitializeEnsembleIndependentMutates() const;

      //! @brief initialize scores that apply to ensembles
      void InitializeFirstEnsembleMutates() const;

      //! @brief initialize scores that apply to two ensembles
      void InitializeTwoEnsembleMutates() const;

      //! @brief initialize scores that will not change for a second ensemble
      void InitializeEnsembleIndependentScores() const;

      //! @brief initialize scores that will not change for a second ensemble
      void InitializeFirstEnsembleScores() const;

      //! @brief initialize scores that depend on two ensembles
      void InitializeTwoEnsembleScores() const;

      //! @brief gives the list of residues which the optimization will be based on
      //! @return reference to list of residues that the optimization will be based on
      const util::ShPtrVector< biol::AABase> &GetAAList() const;

      //! @brief gives chains from fasta command line flag
      //! return shptr list of chains created from fasta command line flag
      const util::ShPtrList< assemble::Chain> &GetChainsFromFastas() const;

      //! @brief gives the set of residue types that should not be considered
      //! @return the set of residue types that should not be considered
      const storage::Set< biol::AAType> &GetExcludedAATypes() const;

      //! @brief gives the sse pool created from the sse pool file given over the command line
      //! @return shptr to sse pool created from the sse pool file given over the command line
      const util::ShPtr< assemble::SSEPool> &GetSSEPool() const;

      //! @brief gives the complete data set as created from the AAList
      //! @return all data pairs created from AAList
      const util::ShPtr< restraint::DataSetPairwise> &GetCompleteDataSet() const;

      //! @brief calculates the distance changes between two ensembles for data pairs given by a data set
      //! @param ENSEMBLE_A the first ensemble used to calculate distance changes
      //! @param ENSEMBLE_B the second ensemble used to calculate distance changes
      //! @param DATA_SET the data pairs for which distance changes will be calculated
      //! @return list of data pairs sorted from largest to smallest distance change
      const util::ShPtr< storage::List< restraint::DataPairwise> > &GetDistanceChangeSortedData
      (
        const assemble::ProteinEnsemble &ENSEMBLE_A, const assemble::ProteinEnsemble &ENSEMBLE_B,
        const restraint::DataSetPairwise &DATA_SET
      ) const;

      //! @brief prints data set to file
      //! @param DATA_SET the data set that will be printed to file
      //! @PARAM ITERATION the iteration number to print
      void Print( const storage::Pair< restraint::DataSetPairwise, double> &RESULT_PAIR, const size_t ITERATION, const std::string &TAG) const;

      //! @brief prints data set to file
      //! @param DATA_SET the data set that will be printed to file
      void PrintDataSet( const storage::Pair< restraint::DataSetPairwise, double> &RESULT_PAIR, const size_t ITERATION, const std::string &TAG) const;

      //! @brief write the data set formatted in to bcl restraints and with distance calculated from given structures
      //! @param DATA_SET the data set for which restraints will be created
      //! @param ITERATION the round of optimization that this data set was made for
      //! @param TAG the tag that should be used for the outputted file
      void WriteRestraints
      (
        const restraint::DataSetPairwise &DATA_SET, const size_t ITERATION, const std::string &TAG
      ) const;

      //! @brief write the data set formatted in to rosetta restraints and with distance calculated from given structures
      //! @param DATA_SET the data set for which restraints will be created
      //! @param ITERATION the round of optimization that this data set was made for
      //! @param TAG the tag that should be used for the outputted file
      void WriteRestraintsRosettaFormat
      (
        const restraint::DataSetPairwise &DATA_SET, const size_t ITERATION, const std::string &TAG
      ) const;

      //! @brief writes pymol script formatted file to show distances on a structure
      void ShowDistancesInPymol
      (
        const storage::Pair< restraint::DataSetPairwise, double> &RESULT_PAIR, const size_t ITERATION, const std::string &TAG
      ) const;

      //! @brief writes pymol script formatted file to show distances on a structure
      void ShowDataSetDistancesInPymol
      (
        const storage::Pair< restraint::DataSetPairwise, double> &RESULT_PAIR, const size_t ITERATION, const std::string &TAG
      ) const;

      //! @brief gives the ensembles that will be used
      //! @return vector of pointers to ensembles
      const util::ShPtrVector< assemble::ProteinEnsemble> &GetEnsembles() const;

      static const ApplicationType OptimizeDataSetPairwise_Instance;

    }; // class OptimizeDataSetPairwise

    //! @brief the Main function
    //! @return error code - 0 for success
    int OptimizeDataSetPairwise::Main() const
    {
      // initialize scores and mutates
      InitializeScores();
      InitializeMutates();

      // true if user wants possible scores and mutates written to files
      if( GetFlagWritePossibleScores()->GetFlag() || GetFlagWritePossibleMutates()->GetFlag())
      {
        const MathFunctionsWeightSet< restraint::DataSetPairwise, double> score_weightset
        (
          GetMathFunctions< restraint::DataSetPairwise, double>(), 0
        );

        const storage::Table< double> score_table( score_weightset.CreateTable());

        if( GetFlagWritePossibleScores()->GetFlag())
        {
          const std::string &filename( GetFlagWritePossibleScores()->GetFirstParameter()->GetValue());
          io::OFStream write;
          io::File::MustOpenOFStream( write, filename);
          score_table.WriteFormatted( write);
        }

        const MathFunctionsWeightSet< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> > mutate_weightset
        (
          GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >(), 0
        );

        const storage::Table< double> mutate_table( mutate_weightset.CreateTable());

        if( GetFlagWritePossibleMutates()->GetFlag())
        {
          const std::string &filename( GetFlagWritePossibleMutates()->GetFirstParameter()->GetValue());
          io::OFStream write;
          io::File::MustOpenOFStream( write, filename);
          mutate_table.WriteFormatted( write);
        }

        // stop if either the mutates or the scores were desired to be printed
        return 0;
      }

      const size_t total_optimizations( fold::DefaultFlags::GetFlagNumberModels()->GetFirstParameter()->GetNumericalValue< size_t>());

      // number of minimizations to be done
      for( size_t current_optimization( 0); current_optimization != total_optimizations; ++current_optimization)
      {
        BCL_MessageCrt
        (
          "optimization #" + util::Format()( current_optimization + 1) + "/" + util::Format()( current_optimization)
        );

        // data set that will be optimized
        util::ShPtr< restraint::DataSetPairwise> data_set( new restraint::DataSetPairwise());

        double data_set_score( 0);

        // mutate the starting data set
        if( GetFlagReadMutatesWeightsStart()->GetFlag())
        {
          const std::string filename( GetFlagReadMutatesWeightsStart()->GetFirstParameter()->GetValue());
          // read in the mutate weightset
          io::IFStream read;
          io::File::MustOpenIFStream( read, filename);
          storage::Table< double> mutate_weight;
          mutate_weight.ReadFormatted( read);
          const MathFunctionsWeightSet< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >
          mutate_weightset( mutate_weight);
          io::File::CloseClearFStream( read);
          util::ShPtr< math::MutateCombine< restraint::DataSetPairwise> > mutate( mutate_weightset.ConstructMutateCombine());
          BCL_MessageStd( "starting Mutate");
          math::MutateResult< restraint::DataSetPairwise> result( mutate->operator()( *data_set));
          BCL_Assert( result.GetArgument().IsDefined(), "result pointer is not defined");
          BCL_MessageStd( "after starting Mutate");
          data_set = result.GetArgument();
          BCL_Assert( data_set.IsDefined(), "data set pointer is not defined");
          BCL_MessageStd( "data_set = result.GetArgument()");
        }

        // true if optimization needs to be done
        if( GetFlagReadScoresWeightsOptimization()->GetFlag() && GetFlagReadScoresWeightsOptimization()->GetFlag())
        {
          BCL_MessageStd( "true if optimization needs to be done");

          // read in the score weightset
          io::IFStream read;
          io::File::MustOpenIFStream( read, GetFlagReadScoresWeightsOptimization()->GetFirstParameter()->GetValue());
          storage::Table< double> score_weights;
          score_weights.ReadFormatted( read);
          const MathFunctionsWeightSet< restraint::DataSetPairwise, double> score_weightset( score_weights);
          io::File::CloseClearFStream( read);
          const util::ShPtr
          <
            math::SumFunction< restraint::DataSetPairwise, double>
          > score( score_weightset.ConstructSumFunction());

          BCL_Assert( score.IsDefined(), "score function for optimization is not defined");

          PrintDataSet
          (
            storage::Pair< restraint::DataSetPairwise, double>( *data_set, score->operator()( *data_set)),
            current_optimization,
            "start_"
          );

          // read in the mutate weightset
          io::File::MustOpenIFStream( read, GetFlagReadMutatesWeightsOptimization()->GetFirstParameter()->GetValue());
          storage::Table< double> mutate_weights;
          mutate_weights.ReadFormatted( read);
          const MathFunctionsWeightSet< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >
            mutate_weightset( mutate_weights);
          const math::ObjectProbabilityDistribution< math::MutateInterface< restraint::DataSetPairwise> > obj_prob_dist
          (
            mutate_weightset.ConstructObjectProbabilityDistribution()
          );
          const math::MutateDecisionNode< restraint::DataSetPairwise> mutate( obj_prob_dist);
          const util::ShPtr< math::MutateInterface< restraint::DataSetPairwise> > sp_mutate( mutate.Clone());

          // create the temperature control
          util::ShPtr< mc::TemperatureInterface> sp_temperature
          (
            new mc::TemperatureAccepted
            (
              // start fraction
              mc::TemperatureAccepted::GetParameterStartFraction()->GetNumericalValue< double>(),
              // end fraction
              mc::TemperatureAccepted::GetParameterEndFraction()->GetNumericalValue< double>(),
              // total number of steps
              fold::DefaultFlags::GetFlagMCNumberIterations()->GetParameterList()( 0)->GetNumericalValue< size_t>(),
              // start temperature
              mc::TemperatureAccepted::GetParameterStartTemperature()->GetNumericalValue< double>(),
              // nr steps between each update
              mc::TemperatureAccepted::GetParameterUpdateInterval()->GetNumericalValue< size_t>()
            )
          );

          // create the metropolis
          util::ShPtr< mc::Metropolis< double> > sp_metropolis
          (
            new mc::Metropolis< double>( sp_temperature, true, 0.0001)
          );

          // create the termination criteria
          opti::CriterionCombine< restraint::DataSetPairwise, double> criterion_combine;
          criterion_combine.InsertCriteria
          (
            opti::CriterionNumberIterations< restraint::DataSetPairwise, double>
            (
              fold::DefaultFlags::GetFlagMCNumberIterations()->GetParameterList()( 0)->GetNumericalValue< size_t>()
            )
          );
          criterion_combine.InsertCriteria
          (
            opti::CriterionUnimproved< restraint::DataSetPairwise, double>
            (
              fold::DefaultFlags::GetFlagMCNumberIterations()->GetParameterList()( 1)->GetNumericalValue< size_t>()
            )
          );

          // create the approximator
          mc::Approximator< restraint::DataSetPairwise, double> approximator
          (
            *score,
            *sp_mutate,
            *sp_metropolis,
            criterion_combine,
            *data_set
          );

          // run the approximator
          approximator.Approximate();

          *data_set = approximator.GetTracker().GetBest()->First();
          data_set_score = approximator.GetTracker().GetBest()->Second();
          BCL_MessageCrt( "scores are: ");

          score->WriteDetailedSchemeAndValues( *data_set, std::cout);
        }

        // mutate the final dataset
        if( GetFlagReadMutatesWeightsEnd()->GetFlag())
        {
          BCL_MessageStd( "mutate the final dataset");
          const std::string filename( GetFlagReadMutatesWeightsEnd()->GetFirstParameter()->GetValue());
          // read in the mutate weightset
          io::IFStream read;
          io::File::MustOpenIFStream( read, filename);
          storage::Table< double> mutate_weight;
          mutate_weight.ReadFormatted( read);
          const MathFunctionsWeightSet< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >
            mutate_weightset( mutate_weight);
          io::File::CloseClearFStream( read);

          util::ShPtr< math::MutateCombine< restraint::DataSetPairwise> > mutate( mutate_weightset.ConstructMutateCombine());
          math::MutateResult< restraint::DataSetPairwise> result( mutate->operator()( *data_set));

          data_set = result.GetArgument();
        }

        BCL_MessageStd( "print the final dataset");
        BCL_Assert( data_set.IsDefined(), "data set pointer is not defined");
        const storage::Pair< restraint::DataSetPairwise, double> result_pair( *data_set, data_set_score);
        if( GetEnsembles().GetSize() == 2)
        {
          Print( result_pair, current_optimization, "distance_change_sorted");
        }
        PrintDataSet( result_pair, current_optimization, "final");

        if( GetFlagRestraintDistanceEnsemble()->GetFlag())
        {
          if( m_WriteRestraintsBCLFormat->GetFlag())
          {
            WriteRestraints( *data_set, current_optimization, "final");
          }
          else
          {
            WriteRestraintsRosettaFormat( *data_set, current_optimization, "final");
          }
        }

        if( GetFlagPymolOutput()->GetFlag())
        {
          //ShowDistancesInPymol( result_pair, 0, "final");
          ShowDataSetDistancesInPymol( result_pair, current_optimization, "final");
        }
      }

      return 0;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &OptimizeDataSetPairwise::GetReadMe() const
    {
      static const std::string s_readme( "readme");
      return s_readme;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> OptimizeDataSetPairwise::InitializeCommand() const
    {
      // initialize ShPtr to a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // flag for fasta
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagFastaRead());

      // flag for minimum sse sizes in protein models
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

      // flag for minimum sse sizes in pool
      sp_cmd->AddFlag( assemble::SSEPool::GetFlagMinSSELengths());

      // flag for reading sse pool
      sp_cmd->AddFlag( assemble::SSEPool::GetFlagPoolRead());

      // aa type exclusion
      sp_cmd->AddFlag( GetFlagAATypeExclusions());

      // desired possible size range of the data set
      sp_cmd->AddFlag( GetFlagDataSetSizeRange());

      // files containing coordinates for exclusion
      sp_cmd->AddFlag( GetFlagCoordinateExclusionFiles());

      // the x,y,z columns containing coordinates for exclusion
      sp_cmd->AddFlag( GetFlagCoordinateExclusionFileColumns());

      // radius around exclusion coordinates that will be used to exclude residues
      sp_cmd->AddFlag( GetFlagCoordinateExclusionRadius());

      // files containing lists of pdbs to make ensembles
      sp_cmd->AddFlag( GetFlagEnsembleFiles());

      // the radius around residues that will keep close-by residues from being considered
      sp_cmd->AddFlag( GetFlagTriangulationRadius());

      // the min and max distances that are desired to be possible for the dataset
      sp_cmd->AddFlag( GetFlagEuclidianDistance());

      // flag to specify desired residue exposure
      sp_cmd->AddFlag( GetFlagExposure());

      // flag to write scores
      sp_cmd->AddFlag( GetFlagWritePossibleScores());

      // flag to write mutates
      sp_cmd->AddFlag( GetFlagWritePossibleMutates());

      // flag to read mutates for starting data set
      sp_cmd->AddFlag( GetFlagReadMutatesWeightsStart());

      // flag to read scores for optimization
      sp_cmd->AddFlag( GetFlagReadScoresWeightsOptimization());

      // flag to read scores for optimization
      sp_cmd->AddFlag( GetFlagReadMutatesWeightsOptimization());

      // flag to read mutates to apply to data set after optimization
      sp_cmd->AddFlag( GetFlagReadMutatesWeightsEnd());

      // flag for monte carlo minimization, max number of rejected steps and max iterations
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagMCNumberIterations());

      // flag for setting temperature adjustment based on accepted steps ratio
      sp_cmd->AddFlag( mc::TemperatureAccepted::GetFlagTemperature());

      // flag for specifying a prefix to be used for writing files
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagPrefix());

      // flag for writing a pymol script to show the distances
      sp_cmd->AddFlag( GetFlagPymolOutput());

      // flag for specifying the number of minimizations that should be done
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagNumberModels());

      // flag for specifying the structures that should be used for calculating distances for restraint output
      sp_cmd->AddFlag( GetFlagRestraintDistanceEnsemble());

      // flag for specifying number of restraints as fraction of number of residue in sses
      sp_cmd->AddFlag( GetFlagDataSetSizeFractionOfPool());

      // flag for specifying number of restraints as connecting all sses a given number of times
      sp_cmd->AddFlag( GetFlagDataSetSizeConnectSSEsInPool());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! default constructor
    OptimizeDataSetPairwise::OptimizeDataSetPairwise() :
      m_WriteRestraintsBCLFormat
      (
        new command::FlagStatic
        (
          "write_bcl_format",
          "If the ensemble of structures flag is provided and this flag is set, restraints will be output in "
          "bcl format. If this flag is not set and an ensemble is provided, the restraints will be output in Rosetta"
          " format."
        )
      )
    {
    }

    //! @brief initializes and fills the enumerator with all scores
    void OptimizeDataSetPairwise::InitializeScores() const
    {
      InitializeEnsembleIndependentScores();
      InitializeFirstEnsembleScores();
      InitializeTwoEnsembleScores();
    }

    //! @brief initializes and fills the enumerator with all mutates
    void OptimizeDataSetPairwise::InitializeMutates() const
    {
      InitializeEnsembleIndependentMutates();
      InitializeFirstEnsembleMutates();
      InitializeTwoEnsembleMutates();
    }

    //! @brief initialize scores that will not change for a second ensemble
    void OptimizeDataSetPairwise::InitializeEnsembleIndependentMutates() const
    {
      // initialize MutateDataSetPairwiseAdd add all
      {
        const Mutate mutate
        (
          new restraint::MutateDataSetPairwiseAdd
          (
            GetCompleteDataSet(), GetCompleteDataSet()->GetSize(), GetCompleteDataSet()->GetSize(),
            restraint::MutateDataSetPairwiseAdd::GetDefaultScheme() + "_all"
          )
        );
        BCL_MessageStd( "To enum adding " + mutate->GetScheme());
        GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( mutate);
      }

      // initialize MutateDataSetPairwiseAdd add one
      {
        const Mutate mutate
        (
          new restraint::MutateDataSetPairwiseAdd
          (
            GetCompleteDataSet(), restraint::MutateDataSetPairwiseAdd::GetDefaultScheme() + "_single"
          )
        );

        BCL_MessageStd( "To enum adding " + mutate->GetScheme());
        GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( mutate);
      }

      // initialize MutateDataSetPairwiseFilterAAType
      {
        const Mutate mutate
        (
          new restraint::MutateDataSetPairwiseFilterAAType( GetExcludedAATypes())
        );

        BCL_MessageStd( "To enum adding " + mutate->GetScheme());
        GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( mutate);
      }

      // initialize MutateDataSetPairwiseFilterSSESize
      if( io::DirectoryEntry( assemble::SSEPool::GetFlagPoolRead()->GetFirstParameter()->GetValue()).DoesExist())
      {
        const Mutate mutate
        (
          new restraint::MutateDataSetPairwiseFilterSSESize
          (
            GetSSEPool(), assemble::SSEPool::GetCommandLineMinSSELengths()
          )
        );

        BCL_MessageStd( "To enum adding " + mutate->GetScheme());
        GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( mutate);
      }

      // initialize MutateDataSetPairwiseRemove
      {
        const Mutate mutate
        (
          new restraint::MutateDataSetPairwiseRemove
          (
            restraint::MutateDataSetPairwiseRemove::GetDefaultScheme() + "_single"
          )
        );

        BCL_MessageStd( "To enum adding " + mutate->GetScheme());
        GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( mutate);
      }

      // initialize mutate for swapping
      {
        // make list to hold the add and remove moves that will make up the swap move
        util::ShPtrList< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> > > mutates;
        mutates.PushBack
        (
          util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> > >
          (
            new restraint::MutateDataSetPairwiseRemove
            (
              "swap_remove"
            )
          )
        );
        mutates.PushBack
        (
          util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> > >
          (
            new restraint::MutateDataSetPairwiseAdd
            (
              GetCompleteDataSet(), "swap_add"
            )
          )
        );

        const Mutate mutate( new math::MutateCombine< restraint::DataSetPairwise>( mutates, false, "swap"));

        BCL_MessageStd( "To enum adding " + mutate->GetScheme());
        GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( mutate);
      }
    }

    //! @brief initialize scores that apply to ensembles
    void OptimizeDataSetPairwise::InitializeFirstEnsembleMutates() const
    {

      const util::ShPtrVector< command::ParameterInterface> &param_list( GetFlagEnsembleFiles()->GetParameterList());

      // initialize mutate for each ensemble given
      for( size_t ensemble_num( 0); ensemble_num < GetEnsembles().GetSize(); ++ensemble_num)
      {
        // initialize DataSetPairwiseCoordinateExclusion
        if( GetFlagCoordinateExclusionFiles()->GetParameterList().GetSize() == GetEnsembles().GetSize())
        {
          io::IFStream read;
          io::File::MustOpenIFStream( read, GetFlagCoordinateExclusionFiles()->GetParameterList()( ensemble_num)->GetValue());
          const storage::Vector< double> columns( GetFlagCoordinateExclusionFileColumns()->GetNumericalList< double>());

          const Mutate score
          (
            new restraint::MutateDataSetPairwiseFilterCoordinateExclusion
            (
              GetFlagCoordinateExclusionRadius()->GetFirstParameter()->GetNumericalValue< double>(),
              read,
              columns( 0), columns( 1), columns( 2),
              GetCompleteDataSet()->GetUniqueDataPoints(),
              *GetEnsembles()( ensemble_num),
              restraint::MutateDataSetPairwiseFilterCoordinateExclusion::GetDefaultScheme() + "_" + util::Format()( ensemble_num)
            )
          );

          GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( score);
        }
        else //< warn user
        {
          BCL_MessageStd
          (
            "Did not initialize restraint::MutateDataSetPairwiseFilterCoordinateExclusion since number of exclusion files " +
            util::Format()( GetFlagCoordinateExclusionFiles()->GetParameterList().GetSize()) +
            " does not match number of ensemble files given " + util::Format()( param_list.GetSize())

          );
        }

        // initialize MutateDataSetPairwiseFilterEuclidianDistance
        {
          const size_t min_param( 0);
          const size_t max_param( 1);
          const Mutate mutate
          (
            new restraint::MutateDataSetPairwiseFilterEuclidianDistance
            (
              math::Range< double>
              (
                GetFlagEuclidianDistance()->GetParameterList()( min_param)->GetNumericalValue< double>(),
                GetFlagEuclidianDistance()->GetParameterList()( max_param)->GetNumericalValue< double>()
              ),
              GetEnsembles()( ensemble_num),
              restraint::MutateDataSetPairwiseFilterEuclidianDistance::GetDefaultScheme() + "_" + util::Format()( ensemble_num)
            )
          );

          GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( mutate);
        }

        // initialize MutateDataSetPairwiseFilterExposure
        {
          const Mutate mutate
          (
            new restraint::MutateDataSetPairwiseFilterExposure
            (
              GetFlagExposure()->GetFirstParameter()->GetNumericalValue< double>(),
              *GetEnsembles()( ensemble_num),
              restraint::MutateDataSetPairwiseFilterExposure::GetDefaultScheme() + "_" + util::Format()( ensemble_num)
            )
          );

          GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( mutate);
        }
      }
    } //< end InitializeFirstEnsembleMutates

    //! @brief initialize scores that apply to two ensembles
    void OptimizeDataSetPairwise::InitializeTwoEnsembleMutates() const
    {
      const util::ShPtrVector< command::ParameterInterface> &param_list( GetFlagEnsembleFiles()->GetParameterList());

      // true if not the right number of ensembles provided
      if( param_list.GetSize() < 2)
      {
        return;
      }

      // initialize MutateDataSetPairwiseFilterTriangulation for first ensemble
      {
        BCL_MessageStd( "initializing MutateDataSetPairwiseFilterTriangulation");
        const assemble::ProteinEnsemble &ensemble_a( *GetEnsembles()( 0));
        const assemble::ProteinEnsemble &ensemble_b( *GetEnsembles()( 1));
        const storage::List< restraint::DataPairwise> sorted_data
        (
          *GetDistanceChangeSortedData( ensemble_a, ensemble_b, *GetCompleteDataSet())
        );
        BCL_MessageStd( "finished GetDistanceChangeSortedData");
        // for first ensemble
        {
          const Mutate mutate
          (
            new restraint::MutateDataSetPairwiseFilterTriangulation
            (
              GetFlagTriangulationRadius()->GetFirstParameter()->GetNumericalValue< double>(),
              GetEnsembles()( 0),
              sorted_data,
              restraint::MutateDataSetPairwiseFilterTriangulation::GetDefaultScheme() + "_" + util::Format()( 0)
            )
          );
          GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( mutate);
        }
        // for second ensemble
        {
          const Mutate mutate
          (
            new restraint::MutateDataSetPairwiseFilterTriangulation
            (
              GetFlagTriangulationRadius()->GetFirstParameter()->GetNumericalValue< double>(),
              GetEnsembles()( 1),
              sorted_data,
              restraint::MutateDataSetPairwiseFilterTriangulation::GetDefaultScheme() + "_" + util::Format()( 1)
            )
          );
          GetMathFunctions< restraint::DataSetPairwise, math::MutateResult< restraint::DataSetPairwise> >().AddFunction( mutate);
        }
      }
    }

    //! @brief initialize scores that will not change for a second ensemble
    void OptimizeDataSetPairwise::InitializeEnsembleIndependentScores() const
    {
       // initialize DataSetPairwiseDataDensity
       {
         const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
         (
           new score::DataSetPairwiseDataDensity( GetAAList().GetSize())
         );

         GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
       }

       // initialize DataSetPairwiseResidueTypeExclusion
       {
         const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
         (
           new score::DataSetPairwiseResidueTypeExclusion( GetExcludedAATypes())
         );

         GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
       }

       // initialize DataSetPairwiseSequenceSeparation
       {
         const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
         (
           new score::DataSetPairwiseSequenceSeparation( GetAAList().GetSize())
         );

         GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
       }

       // initialize DataSetPairwiseSize
       {
         // get min and max size as specified by flag
         size_t min_size( GetFlagDataSetSizeRange()->GetNumericalList< size_t>()( 0));
         size_t max_size( GetFlagDataSetSizeRange()->GetNumericalList< size_t>()( 1));

         // if the user wants size of dataset to be a fraction of residues in the pool
         if( GetFlagDataSetSizeFractionOfPool()->GetFlag())
         {
           // fraction of residues in sse pool desired by user
           const double fraction( GetFlagDataSetSizeFractionOfPool()->GetFirstParameter()->GetNumericalValue< double>());

           // make sure the pool file exists
           if( io::DirectoryEntry( assemble::SSEPool::GetFlagPoolRead()->GetFirstParameter()->GetValue()).DoesExist())
           {
             // get random non over lapping set of sses from the pool
             const storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> sses
             (
               GetSSEPool()->GetRandomNonOverlappingSet()
             );

             // to hold the number of residues in the sses
             size_t num_sse_residues( 0);

             // iterate over the sses to sum up the number of residues
             for
             (
               storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
                 sse_itr( sses.Begin()), sse_itr_end( sses.End());
               sse_itr != sse_itr_end;
               ++sse_itr
             )
             {
               num_sse_residues += ( *sse_itr)->GetSize();
             }

             BCL_MessageStd
             (
               "number of residues is " + util::Format()( num_sse_residues) +
               " fraction desired of restraints is " + util::Format()( fraction)
             );

             // set min and max data set size
             min_size = double( num_sse_residues) * fraction;
             max_size = double( num_sse_residues) * fraction;
           }
           else
           {
             BCL_MessageCrt
             (
               "pool file doesn't exist. can't set number of restraint as fraction of residues in sses"
             );
           }
         }
         // if the user wants size of dataset to be size of connecting each SSE a given number of times
         if( GetFlagDataSetSizeConnectSSEsInPool()->GetFlag())
         {
           // the number of times each sse should be connected
           const size_t connections
           (
             GetFlagDataSetSizeConnectSSEsInPool()->GetFirstParameter()->GetNumericalValue< size_t>()
           );

           // make sure the pool file exists
           if( io::DirectoryEntry( assemble::SSEPool::GetFlagPoolRead()->GetFirstParameter()->GetValue()).DoesExist())
           {
             // get random non over lapping set of sses from the pool
             const storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> sses
             (
               GetSSEPool()->GetRandomNonOverlappingSet()
             );

             // to hold the number of residues in the sses
             const size_t num_sses( sses.GetSize());

             // number of needed restraints
             const size_t num_restraints( connections * num_sses * ( num_sses - 1) / 2);

             BCL_MessageStd
             (
               "number of sses is " + util::Format()( num_sses) +
               " and number of connections desired is " + util::Format()( connections) +
               " so number of restraints is therefore " + util::Format()( num_restraints)
             );

             // set min and max data set size
             min_size = num_restraints;
             max_size = num_restraints;
           }
           else
           {
             BCL_MessageCrt
             (
               "pool file doesn't exist. can't set number of restraints to connect sses some number of times"
             );
           }
         }

         // command line cutoffs for min and max override other options
         const size_t min_size_cutoff( GetFlagDataSetSizeRange()->GetNumericalList< size_t>()( 0));
         const size_t max_size_cutoff( GetFlagDataSetSizeRange()->GetNumericalList< size_t>()( 1));

         // set min size and max size if they are outside the range allowed by min_size_cutoff and max_size_cutoff
         if( min_size < min_size_cutoff)
         {
           min_size = min_size_cutoff;
         }
         if( max_size > max_size_cutoff)
         {
           max_size = max_size_cutoff;
         }
         // make sure min and max are meaningfully set
         if( min_size > max_size)
         {
           min_size = max_size;
         }
         if( max_size < min_size)
         {
           max_size = min_size;
         }

         BCL_MessageStd
         (
           "desired data set size is between " + util::Format()( min_size) + " and " +
           util::Format()( max_size)
         );
         const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
         (
           new score::DataSetPairwiseSize
           (
             min_size,
             max_size
           )
         );

         GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
       }

       // initialize DataSetPairwiseSSEConnection
       if( io::DirectoryEntry( assemble::SSEPool::GetFlagPoolRead()->GetFirstParameter()->GetValue()).DoesExist())
       {
         {
           const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
           (
             new score::DataSetPairwiseSSEConnection
             (
               GetSSEPool()
             )
           );

           GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
         }

         // initialize DataSetPairwiseSSESize
         {
           const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
           (
             new score::DataSetPairwiseSSESize
             (
               GetSSEPool()
             )
           );

           GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
         }

         // initialize DataSetPairwiseSSETerm
         {
           const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
           (
             new score::DataSetPairwiseSSETerm
             (
               GetSSEPool()
             )
           );

           GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
         }

         // initialize DataSetPairwiseBipolar
         {
           const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
           (
             new score::DataSetPairwiseBipolar( *GetSSEPool())
           );

           GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
         }

         // initialize DataSetPairwiseSSECenter
         {
           const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
           (
             new score::DataSetPairwiseSSECenter( *GetSSEPool())
           );

           GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
         }
       }
    }

    //! @brief initialize scores that use one ensemble
    void OptimizeDataSetPairwise::InitializeFirstEnsembleScores() const
    {
      // initialize score for each ensemble given
      for( size_t ensemble_num( 0); ensemble_num < GetEnsembles().GetSize(); ++ensemble_num)
      {
        // initialize DataSetPairwiseCoordinateExclusion
        if( GetFlagCoordinateExclusionFiles()->GetParameterList().GetSize() == GetEnsembles().GetSize())
        {
          io::IFStream read;
          io::File::MustOpenIFStream( read, GetFlagCoordinateExclusionFiles()->GetParameterList()( ensemble_num)->GetValue());
          const storage::Vector< double> columns( GetFlagCoordinateExclusionFileColumns()->GetNumericalList< double>());

          const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
          (
            new score::DataSetPairwiseCoordinateExclusion
            (
              GetFlagCoordinateExclusionRadius()->GetFirstParameter()->GetNumericalValue< double>(),
              read,
              columns( 0), columns( 1), columns( 2),
              GetCompleteDataSet()->GetUniqueDataPoints(),
              *GetEnsembles()( ensemble_num),
              score::DataSetPairwiseCoordinateExclusion::GetDefaultScheme() + "_" + util::Format()( ensemble_num)
            )
          );

          GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
        }
        else //< warn user
        {
          BCL_MessageStd
          (
            "Did not initialize DataSetPairwiseCoordinateExclusion since number of exclusion files " +
            util::Format()( GetFlagCoordinateExclusionFiles()->GetParameterList().GetSize()) +
            " does not match number of ensemble files given " + util::Format()( GetEnsembles().GetSize())

          );
        }

        // initialize DataSetPairwiseCoordinateTriangulation
        {
          BCL_MessageStd( "initializing DataSetPairwiseCoordinateTriangulation");
          const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
          (
            new score::DataSetPairwiseCoordinateTriangulation
            (
              GetFlagTriangulationRadius()->GetFirstParameter()->GetNumericalValue< double>(),
              GetEnsembles()( ensemble_num),
              score::DataSetPairwiseCoordinateTriangulation::GetDefaultScheme() + "_" + util::Format()( ensemble_num)
            )
          );

          GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
        }

        // initialize DataSetPairwiseEuclidianDistance
        {
          const size_t min_param( 0);
          const size_t max_param( 1);
          const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
          (
            new score::DataSetPairwiseEuclidianDistance
            (
              math::Range< double>
              (
                GetFlagEuclidianDistance()->GetParameterList()( min_param)->GetNumericalValue< double>(),
                GetFlagEuclidianDistance()->GetParameterList()( max_param)->GetNumericalValue< double>()
              ),
              GetEnsembles()( ensemble_num),
              score::DataSetPairwiseEuclidianDistance::GetDefaultScheme() + "_" + util::Format()( ensemble_num)
            )
          );

          GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
        }

        // initialize DataSetPairwiseStructuralExposure
        {
          const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
          (
            new score::DataSetPairwiseStructuralExposure
            (
              GetFlagExposure()->GetFirstParameter()->GetNumericalValue< double>(),
              *GetEnsembles()( ensemble_num),
              score::DataSetPairwiseStructuralExposure::GetDefaultScheme() + "_" + util::Format()( ensemble_num)
            )
          );

          GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
        }
      }
    } //< end InitializeFirstEnsembleScores

    //! @brief initialize scores that depend on two ensembles
    void OptimizeDataSetPairwise::InitializeTwoEnsembleScores() const
    {
      const util::ShPtrVector< command::ParameterInterface> &param_list( GetFlagEnsembleFiles()->GetParameterList());

      // true if not the right number of ensembles provided
      if( param_list.GetSize() < 2)
      {
        return;
      }

      // initialize DataSetPairwiseDistanceChangeMagnitude
      {
        const util::ShPtr< math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double> > score
        (
          new score::DataSetPairwiseDistanceChangeMagnitude
          (
            *GetEnsembles()( 0),
            *GetEnsembles()( 1),
            *GetCompleteDataSet()
          )
        );

        GetMathFunctions< restraint::DataSetPairwise, double>().AddFunction( score);
      }
    }

    //! @brief gives the list of residues which the optimization will be based on
    //! @return reference to list of residues that the optimization will be based on
    const util::ShPtrVector< biol::AABase> &OptimizeDataSetPairwise::GetAAList() const
    {
      static util::ShPtrVector< biol::AABase> s_residues;

      if( s_residues.IsEmpty())
      {
        for
        (
          util::ShPtrList< assemble::Chain>::const_iterator
            chain_itr( GetChainsFromFastas().Begin()), chain_itr_end( GetChainsFromFastas().End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          s_residues.Append( ( *chain_itr)->GetSequence()->GetData());
        }
      }

      return s_residues;
    }

    //! @brief gives chains from fasta command line flag
    //! return shptr list of chains created from fasta command line flag
    const util::ShPtrList< assemble::Chain> &OptimizeDataSetPairwise::GetChainsFromFastas() const
    {
      static util::ShPtrList< assemble::Chain> s_chains;

      if( s_chains.IsEmpty())
      {
        // iterate through the list of fastas
        for
        (
          util::ShPtrVector< command::ParameterInterface>::const_iterator
            fasta_itr( fold::DefaultFlags::GetFlagFastaRead()->GetParameterList().Begin()),
            fasta_itr_end( fold::DefaultFlags::GetFlagFastaRead()->GetParameterList().End());
          fasta_itr != fasta_itr_end;
          ++fasta_itr
        )
        {
          // remove fasta extension
          const std::string fasta_id( io::File::RemoveFullExtension( ( *fasta_itr)->GetValue()));

          // get the chain id, which should be the last char before the ".fasta" extension
          const char chain_id( fasta_id[ fasta_id.length() - 1]);
          // open fasta file
          io::IFStream read;
          io::File::MustOpenIFStream( read, ( *fasta_itr)->GetValue());
          // insert chain for fasta file
          s_chains.PushBack
          (
            pdb::Factory( biol::GetAAClasses().e_AABackBone).ChainFromFastaStream( chain_id, read)
          );
        }
      }

      return s_chains;
    }

    //! @brief gives the set of residue types that should not be considered
    //! @return the set of residue types that should not be considered
    const storage::Set< biol::AAType> &OptimizeDataSetPairwise::GetExcludedAATypes() const
    {
      static const storage::Set< biol::AAType> s_excluded_aa_types
      (
        GetFlagAATypeExclusions()->GetObjectSet< biol::AAType>()
      );

      return s_excluded_aa_types;
    }

    //! @brief gives the sse pool created from the sse pool file given over the command line
    //! @return shptr to sse pool created from the sse pool file given over the command line
    const util::ShPtr< assemble::SSEPool> &OptimizeDataSetPairwise::GetSSEPool() const
    {
      static util::ShPtr< assemble::SSEPool> s_pool;

      if( !s_pool.IsDefined())
      {
        s_pool = util::ShPtr< assemble::SSEPool>( new assemble::SSEPool());

        assemble::ProteinModel empty_model;

        // initialize map to hold the min pool lengths
        storage::Map< biol::SSType, size_t> min_pool_sse_lengths( assemble::SSEPool::GetCommandLineMinSSELengths());

        // iterate through the chains to put them in the model
        for
        (
          util::ShPtrList< assemble::Chain>::const_iterator
            chain_itr( GetChainsFromFastas().Begin()), chain_itr_end( GetChainsFromFastas().End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          empty_model.Insert( *chain_itr);
        }

        // read the pool
        // open pool file
        io::IFStream read;
        io::File::MustOpenIFStream( read, assemble::SSEPool::GetFlagPoolRead()->GetFirstParameter()->GetValue());

        // read pool
        s_pool->ReadSSEPool
        (
          read,
          empty_model,
          min_pool_sse_lengths[ biol::GetSSTypes().HELIX],
          min_pool_sse_lengths[ biol::GetSSTypes().STRAND]
        );

        // close stream
        io::File::CloseClearFStream( read);
      }

      return s_pool;
    }

    //! @brief gives the complete data set as created from the AAList
    //! @return all data pairs created from AAList
    const util::ShPtr< restraint::DataSetPairwise> &OptimizeDataSetPairwise::GetCompleteDataSet() const
    {
      static const util::ShPtr< restraint::DataSetPairwise> s_dataset( restraint::DataSetPairwise::GetCompleteDataSet( GetAAList()));

      return s_dataset;
    }

    //! @brief calculates the distance changes between two ensembles for data pairs given by a data set
    //! @return list of data pairs sorted from largest to smallest magnitude distance change
    const util::ShPtr< storage::List< restraint::DataPairwise> > &OptimizeDataSetPairwise::GetDistanceChangeSortedData
    (
      const assemble::ProteinEnsemble &ENSEMBLE_A, const assemble::ProteinEnsemble &ENSEMBLE_B,
      const restraint::DataSetPairwise &DATA_SET
    ) const
    {
      static util::ShPtr< storage::List< restraint::DataPairwise> > s_largest_to_smallest_sorted_values;

      if( !s_largest_to_smallest_sorted_values.IsDefined())
      {
        // to hold data pairs sorted from largest to smallest distance change;
        storage::List< restraint::DataPairwise> sorted_values;

                // get the distance changes
        const std::multimap
        <
          math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean
        > distance_changes( ENSEMBLE_A.GetDistanceChangesMeanSD( DATA_SET, ENSEMBLE_B));

        // print the distance changes sorted by magnitude change
        for
        (
          std::multimap< math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean>::const_iterator
            itr( distance_changes.begin()), itr_end( distance_changes.end());
          itr != itr_end;
          ++itr
        )
        {
          BCL_MessageDbg
          (
            "distance change for " + itr->second.GetIdentification() + " is " +
            util::Format()( itr->first.GetAverage())
          );
        }

        // reverse iterate through distance changes to fill sorted_values with distance changes from largest to smallest
        for
        (
          std::multimap
          <
            math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean
          >::const_reverse_iterator
          itr( distance_changes.rbegin()),
          itr_end( distance_changes.rend());
          itr != itr_end; ++itr
        )
        {
          if( itr->first.GetWeight())
          {
            sorted_values.PushBack( itr->second);
          }
        }

        s_largest_to_smallest_sorted_values = util::ShPtr< storage::List< restraint::DataPairwise> >( sorted_values.Clone());
      }

      return s_largest_to_smallest_sorted_values;
    }

    //! @brief return command line flag for specifying the residues that should not be considered
    //! @return command line flag for specifying the residues that should not be considered
    util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagAATypeExclusions()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "exclude_residue_types", "list of residue types that should be excluded",
           command::Parameter
          (
            "residue_three_letter_code",
            "any number of three letter codes from the list to be excluded",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>( biol::GetAATypes().Begin(), biol::GetAATypes().ASX.GetIterator())
            )
          ), 0, biol::AATypes::s_NumberStandardAATypes
        )
      );

      return s_flag;
    }

    //! @brief prints data set to file
    //! @param DATA_SET the data set that will be printed to file
    void OptimizeDataSetPairwise::Print( const storage::Pair< restraint::DataSetPairwise, double> &RESULT_PAIR, const size_t ITERATION, const std::string &TAG) const
    {
      const std::string prefix( fold::DefaultFlags::GetFlagPrefix()->GetFirstParameter()->GetValue());
      const std::string filename
      (
        prefix + "_" + TAG + "_" + mc::PrintInterface< restraint::DataSetPairwise, double>::GetRoundNumberFormat()( ITERATION) + ".data"
      );
      io::OFStream write;

      io::File::MustOpenOFStream( write, filename);

      // print score
      write << RESULT_PAIR.Second() << '\n';

      // get the data set sorted by distance change magnitude
      const std::multimap< math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean>
        magnitude_sorted( GetEnsembles()( 0)->GetDistanceChangesMeanSD( RESULT_PAIR.First(), *GetEnsembles()( 1)));

      // iterate through the sorted data
      for
      (
        std::multimap< math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean>::const_iterator
        itr( magnitude_sorted.begin()), itr_end( magnitude_sorted.end());
        itr != itr_end; ++itr
      )
      {
        write << itr->second.GetIdentification() << " mean_dist_change " << itr->first.GetAverage() << " std_dev "
              << itr->first.GetStandardDeviation() << '\n';
      }
      io::File::CloseClearFStream( write);
    }

    //! @brief prints data set to file
    //! @param DATA_SET the data set that will be printed to file
    void OptimizeDataSetPairwise::PrintDataSet
    (
      const storage::Pair< restraint::DataSetPairwise, double> &RESULT_PAIR, const size_t ITERATION, const std::string &TAG
    ) const
    {
      const std::string prefix( fold::DefaultFlags::GetFlagPrefix()->GetFirstParameter()->GetValue());
      const std::string filename
      (
        prefix + "_" + TAG + "_" + mc::PrintInterface< restraint::DataSetPairwise, double>::GetRoundNumberFormat()( ITERATION) + ".data"
      );
      io::OFStream write;

      io::File::MustOpenOFStream( write, filename);

      // print score
      write << "score: " << RESULT_PAIR.Second() << '\n';

      // iterate through the sorted data
      for
      (
        restraint::DataSetPairwise::const_iterator
        itr( RESULT_PAIR.First().Begin()), itr_end( RESULT_PAIR.First().End());
        itr != itr_end; ++itr
      )
      {
        write << itr->GetIdentification() << '\n';
      }
      io::File::CloseClearFStream( write);
    }

    //! @brief write the data set formatted in to bcl restraints and with distance calculated from given structures
    //! @param DATA_SET the data set for which restraints will be created
    //! @param ITERATION the round of optimization that this data set was made for
    //! @param TAG the tag that should be used for the outputted file
    void OptimizeDataSetPairwise::WriteRestraints
    (
      const restraint::DataSetPairwise &DATA_SET, const size_t ITERATION, const std::string &TAG
    ) const
    {
      // statically initialize the ensemble from the provided file
      static const assemble::ProteinEnsemble ensemble
      (
        GetFlagRestraintDistanceEnsemble()->GetFirstParameter()->GetValue(), 0, biol::GetAAClasses().e_AAComplete
      );

      // create the restraints from the ensemble and given data set
      const util::ShPtrVector< restraint::AtomDistance> restraints
      (
        restraint::HandlerAtomDistanceAssigned::CreateRestraints( ensemble, DATA_SET)
      );

      // create filename and open it
      const std::string prefix( fold::DefaultFlags::GetFlagPrefix()->GetFirstParameter()->GetValue());
      const std::string filename
      (
        prefix + "_" + TAG + "_" +
        mc::PrintInterface< restraint::DataSetPairwise, double>::GetRoundNumberFormat()( ITERATION) + ".bcl_cst"
      );
      BCL_MessageDbg( "writing bcl restraints to " + filename);
      io::OFStream write;
      io::File::MustOpenOFStream( write, filename);
      restraint::HandlerAtomDistanceAssigned::WriteRestraints( write, restraints);
    }

    //! @brief write the data set formatted in to rosetta restraints and with distance calculated from given structures
    //! @param DATA_SET the data set for which restraints will be created
    //! @param ITERATION the round of optimization that this data set was made for
    //! @param TAG the tag that should be used for the outputted file
    void OptimizeDataSetPairwise::WriteRestraintsRosettaFormat
    (
      const restraint::DataSetPairwise &DATA_SET, const size_t ITERATION, const std::string &TAG
    ) const
    {
      // statically initialize the ensemble from the provided file
      static const assemble::ProteinEnsemble ensemble
      (
        GetFlagRestraintDistanceEnsemble()->GetFirstParameter()->GetValue(), 0, biol::GetAAClasses().e_AAComplete
      );

      // create the restraints from the ensemble and given data set
      const util::ShPtrVector< restraint::AtomDistance> restraints
      (
        restraint::HandlerAtomDistanceAssigned::CreateRestraints( ensemble, DATA_SET)
      );

      // create filename and open it
      const std::string prefix( fold::DefaultFlags::GetFlagPrefix()->GetFirstParameter()->GetValue());
      const std::string filename
      (
        prefix + "_" + TAG + "_" +
        mc::PrintInterface< restraint::DataSetPairwise, double>::GetRoundNumberFormat()( ITERATION) + ".rosetta_cst"
      );
      BCL_MessageDbg( "writing bcl restraints to " + filename);
      io::OFStream write;
      io::File::MustOpenOFStream( write, filename);
      restraint::HandlerAtomDistanceAssigned::WriteDistanceRestraintsRosettaFormat( write, restraints);
    }

    void OptimizeDataSetPairwise::ShowDistancesInPymol
    (
      const storage::Pair< restraint::DataSetPairwise, double> &RESULT_PAIR, const size_t ITERATION, const std::string &TAG
    ) const
    {
      // data pair, mean distance, color
      storage::List< storage::Triplet< restraint::DataPairwise, double, linal::Vector3D> > all_distance;

      // get the data set sorted by distance change magnitude
      const std::multimap< math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean>
        magnitude_sorted( GetEnsembles()( 0)->GetDistanceChangesMeanSD( RESULT_PAIR.First(), *GetEnsembles()( 1)));

      // iterate through the sorted data
      for
      (
        std::multimap< math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean>::const_iterator
        itr( magnitude_sorted.begin()), itr_end( magnitude_sorted.end());
        itr != itr_end; ++itr
      )
      {
        const double mean( itr->first.GetAverage());
        util::Color color( util::GetColors().e_Yellow);
        if( mean < 0)
        {
          color = util::GetColors().e_Green;
        }

        const storage::Triplet< restraint::DataPairwise, double, linal::Vector3D> triplet( itr->second, mean, *color);
        all_distance.PushBack( triplet);
      }
      const std::string prefix( fold::DefaultFlags::GetFlagPrefix()->GetFirstParameter()->GetValue());
      const std::string filename
      (
        prefix + "_" + TAG + "_" +
        mc::PrintInterface< restraint::DataSetPairwise, double>::GetRoundNumberFormat()( ITERATION) + ".pml"
      );
      BCL_MessageStd( "writing pymol script to " + filename);
      io::OFStream write;
      io::File::MustOpenOFStream( write, filename);
      BCL_MessageStd( "all_distance size is " + util::Format()( all_distance.GetSize()));
      restraint::ShowDistancesInPymol( write, all_distance);

      io::File::CloseClearFStream( write);
    }

    //! @brief writes pymol script formatted file to show distances on a structure
    void OptimizeDataSetPairwise::ShowDataSetDistancesInPymol
    (
      const storage::Pair< restraint::DataSetPairwise, double> &RESULT_PAIR, const size_t ITERATION, const std::string &TAG
    ) const
    {
      // data pair, mean distance, color
      storage::List< storage::Triplet< restraint::DataPairwise, double, linal::Vector3D> > all_distance;

      // iterate through the sorted data
      for
      (
        restraint::DataSetPairwise::const_iterator
        itr( RESULT_PAIR.First().Begin()), itr_end( RESULT_PAIR.First().End());
        itr != itr_end; ++itr
      )
      {
        const double mean( util::GetUndefinedDouble());
        util::Color color( util::GetColors().e_Yellow);
        if( mean < 0 && util::IsDefined( mean))
        {
          color = util::GetColors().e_Green;
        }

        const storage::Triplet< restraint::DataPairwise, double, linal::Vector3D> triplet( *itr, mean, *color);
        all_distance.PushBack( triplet);
      }

      const std::string prefix( fold::DefaultFlags::GetFlagPrefix()->GetFirstParameter()->GetValue());
      const std::string filename
      (
        prefix + "_" + TAG + "_" +
        mc::PrintInterface< restraint::DataSetPairwise, double>::GetRoundNumberFormat()( ITERATION) + ".pml"
      );
      BCL_MessageStd( "writing pymol script to " + filename);
      io::OFStream write;
      io::File::MustOpenOFStream( write, filename);
      BCL_MessageStd( "all_distance size is " + util::Format()( all_distance.GetSize()));
      restraint::ShowDistancesInPymol( write, all_distance);

      io::File::CloseClearFStream( write);
    }

    //! @brief gives the ensembles that will be used
    //! @return vector of pointers to ensembles
    const util::ShPtrVector< assemble::ProteinEnsemble> &OptimizeDataSetPairwise::GetEnsembles() const
    {
      static util::ShPtrVector< assemble::ProteinEnsemble> s_ensembles;

      if( s_ensembles.IsEmpty())
      {
        // get the parameter list of ensemble filenames
        const util::ShPtrVector< command::ParameterInterface> &param_list( GetFlagEnsembleFiles()->GetParameterList());

        // iterate through the parameter list
        for
        (
          util::ShPtrVector< command::ParameterInterface>::const_iterator
            param_itr( param_list.Begin()), param_itr_end( param_list.End());
          param_itr != param_itr_end;
          ++param_itr
        )
        {
          // create ensmble
          const util::ShPtr< assemble::ProteinEnsemble> ensemble( new assemble::ProteinEnsemble( ( *param_itr)->GetValue(), 0, biol::GetAAClasses().e_AAComplete));

          // add ensemble to the vector
          s_ensembles.PushBack( ensemble);
        }
      }

      return s_ensembles;
    }

    //! @brief return command line flag for specifying the desired possible size range of the data set
    //! @return command line flag for specifying the desired possible size range of the data set
    util::ShPtr< command::FlagStatic> &OptimizeDataSetPairwise::GetFlagDataSetSizeRange()
    {
      // initialize static flag
      static util::ShPtr< command::FlagStatic> s_flag;

      if( !s_flag.IsDefined())
      {
        s_flag = util::ShPtr< command::FlagStatic>
        (
          new command::FlagStatic
          (
            "data_set_size_range", "desired possible range for the data set inclusive [min,max]"
          )
        );

        // parameter for min
        {
          util::ShPtr< command::ParameterInterface> s_param
          (
            new command::Parameter
            (
              "min_size",
              "inclusive minimum desired size for the data set",
              util::Format()( 10)
            )
          );

          s_flag->PushBack( s_param);
        }

        // parameter for max
        {
          util::ShPtr< command::ParameterInterface> s_param
          (
            new command::Parameter
            (
              "max_size",
              "inclusive maximum desired size for the data set",
              util::Format()( 20)
            )
          );

          s_flag->PushBack( s_param);
        }
      }

      return s_flag;
    }

    //! @brief return command line flag for specifying files containing coordinates for exclusion
    //! @return command line flag for specifying coordinates for exclusion
    util::ShPtr< command::FlagDynamic> &OptimizeDataSetPairwise::GetFlagCoordinateExclusionFiles()
    {
      // initialize static flag
      static util::ShPtr< command::FlagDynamic> s_flag
      (
        new command::FlagDynamic
        (
          "coordinate_exclusion",
          "\t filenames containing lists of coordinates which will be used to exclude coordinates in an ensemble. This order must match ensemble file order.",
          command::Parameter
          (
            "coordinates_filename",
            "\tfilename for inputting a list coordinates that will be used to exclude positions in an ensemble",
            command::ParameterCheckFileExistence()
          ),
          0, 2
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for specifying the x,y,z columns containing coordinates for exclusion
    //! @return command line flag for specifying the x,y,z columns with coordinates for exclusion
    util::ShPtr< command::FlagStatic> &OptimizeDataSetPairwise::GetFlagCoordinateExclusionFileColumns()
    {
      // initialize static flag
      static util::ShPtr< command::FlagStatic> s_flag;

      if( !s_flag.IsDefined())
      {
        s_flag = util::ShPtr< command::FlagStatic>
        (
          new command::FlagStatic
          (
            "coordinate_exclusion_columns", "specify columns where the x,y, and z coordinates are found. columns start with 0."
          )
        );

        // parameter for x column
        {
          util::ShPtr< command::ParameterInterface> s_param
          (
            new command::Parameter
            (
              "x_column",
              "column where x coordinates are",
              util::Format()( 0)
            )
          );

          s_flag->PushBack( s_param);
        }

        // parameter for y column
        {
          util::ShPtr< command::ParameterInterface> s_param
          (
            new command::Parameter
            (
              "y_column",
              "column where y coordinates are",
              util::Format()( 1)
            )
          );

          s_flag->PushBack( s_param);
        }

        // parameter for z column
        {
          util::ShPtr< command::ParameterInterface> s_param
          (
            new command::Parameter
            (
              "z_column",
              "column where z coordinates are",
              util::Format()( 2)
            )
          );

          s_flag->PushBack( s_param);
        }
      }

      return s_flag;
    }

    //! @brief specify the radius around exclusion coordinates that will be used to exclude residues
    //! @return command line flag for specifying the radius around exclusion coordinates that will be used
    util::ShPtr< command::FlagStatic> &OptimizeDataSetPairwise::GetFlagCoordinateExclusionRadius()
    {
      // initialize static flag
      static util::ShPtr< command::FlagStatic> s_flag
      (
        new command::FlagStatic
        (
          "coordinate_exclusion_radius",
          "\tany residues within this radius of one of the exclusion coordinates will not be considered",
          command::Parameter( "radius", "\texclusion radius", "10")
        )
      );

      // end
      return s_flag;
    }

    //! @brief return command line flag for specifying files containing lists of pdbs to make ensembles
    //! @return command line flag for specifying files containing lists of pdbs to make ensembles
    util::ShPtr< command::FlagDynamic> &OptimizeDataSetPairwise::GetFlagEnsembleFiles()
    {
      // initialize static flag
      static util::ShPtr< command::FlagDynamic> s_flag
      (
        new command::FlagDynamic
        (
          "ensembles",
          "\t filenames containing lists of pdbs from which ensembles can be created",
          command::Parameter
          (
            "ensemble_filename",
            "\tfilename for inputting a list of pdb filenames",
            command::ParameterCheckFileExistence()
          ),
          0, 2
        )
      );
      // end
      return s_flag;
    }

    //! @brief specify the radius around residues that will keep close-by residues from being considered
    //! @return command line flag for the radius around residues to keep close-by residues from being considered
    util::ShPtr< command::FlagStatic> &OptimizeDataSetPairwise::GetFlagTriangulationRadius()
    {
      // initialize static flag
      static util::ShPtr< command::FlagStatic> s_flag
      (
        new command::FlagStatic
        (
          "triangulation_radius",
          "\tany residues within this radius of one of a previously used residues will not be considered",
          command::Parameter( "radius", "\radius", "15")
        )
      );

      // end
      return s_flag;
    }

    //! @brief specify the min and max distances that are desired to be possible for the dataset
    //! @return command line flag for specifying the min and max distances
    util::ShPtr< command::FlagStatic> &OptimizeDataSetPairwise::GetFlagEuclidianDistance()
    {
      // initialize static flag
      static util::ShPtr< command::FlagStatic> s_flag;

      if( !s_flag.IsDefined())
      {
        s_flag = util::ShPtr< command::FlagStatic>
        (
          new command::FlagStatic
          (
            "distance_min_max", "range of distance measures that are desired"
          )
        );

        // parameter for min
        {
          util::ShPtr< command::ParameterInterface> s_param
          (
            new command::Parameter
            (
              "min",
              "min distance that should be in data set",
              util::Format()( 15)
            )
          );

          s_flag->PushBack( s_param);
        }

        // parameter for max
        {
          util::ShPtr< command::ParameterInterface> s_param
          (
            new command::Parameter
            (
              "max",
              "max distance that should be in data set",
              util::Format()( 50)
            )
          );

          s_flag->PushBack( s_param);
        }
      }

      return s_flag;
    }

    //! @brief return command line flag for specifying desired residue exposure
    //! @return command line flag for specifying desired residue exposure
    util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagExposure()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "nc_limit", "\t\tmaximum neighbor count a residue should have",
          command::Parameter( "neighbor_count", "\tmaximum neighbor count", "8")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for writing scores to a file
    //! @return command line flag for writing scores to a file
    util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagWritePossibleScores()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "write_scores", "\t\twrite the scores that can be possibly used based on given command line options in table format with default weights",
          command::Parameter( "filename", "\tfilename for writing table", "scores_weights.table")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for writing mutates to a file
    //! @return command line flag for writing mutates to a file
    util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagWritePossibleMutates()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "write_mutates", "\t\twrite the mutates that can be possibly used based on given command line options in table format with default weights",
          command::Parameter( "filename", "\tfilename for writing table", "mutates_weights.table")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for reading mutates from a table file for creating starting dataset
    //! @return command line flag for reading mutates from a table file for creating starting dataset
    util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagReadMutatesWeightsStart()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "read_mutates_start", "\t\tread the mutates with weights from table formatted file that will create a data set without optimization. If optimization flag is used, this will be the starting dataset for optimization. Here, the weights indicate the order in which the mutate is applied",
          command::Parameter( "filename", "\tfilename for reading table", "start_mutates_weights.table")
        )
      );
      // end
      return s_flag;
    }

    //! @brief command line flag for reading scores from a table file to use for optimizing the starting data set
    //! @return command line flag for reading mutates from a table file to use for optimizing the starting data set
    util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagReadScoresWeightsOptimization()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "read_scores_optimization", "\t\tread the scores with weights from table formatted file that will be used during optimization",
          command::Parameter( "filename", "\tfilename for writing table", "optimization_scores_weights.table")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for reading mutates from a table file for optimizing the starting data set
    //! @return command line flag for reading mutates from a table file for optimizing the starting data set
    util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagReadMutatesWeightsOptimization()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "read_mutates_optimization",
          "\t\tread the mutates with weights from table formatted file that will be used during optimization. Here, the weights the probability with which a mutate is applied",
          command::Parameter( "filename", "\tfilename for reading table", "optimization_mutates_weights.table")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for reading mutates from a table file for filtering optimized data set
    //! @return command line flag for reading mutates from a table file for filtering optimized data set
    util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagReadMutatesWeightsEnd()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "read_mutates_end",
          "\t\tread the mutates with weights from table formatted file that will be applied to data set after optimization. Here, the weights indicate the order in which the mutate is applied",
          command::Parameter( "filename", "\tfilename for reading table", "end_mutates_weights.table")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for writing a pymol script to show the distances
    //! @return command line flag for writing a pymol script to show the distances
    util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagPymolOutput()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "pymol_output",
          "\t\twrite pymol script to show distances."
        )
      );
      // end
      return s_flag;
    }

    //! @brief the filename containing pdbs of the structures that should be used to calculate the distance restraints
    //! @return command line flag for specifying pdbs to be used for calculated restraint distances
    util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagRestraintDistanceEnsemble()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "restaint_distance_structures",
          "\t\tspecify the filename containing pdbs of the structures that should be used to calculate the actual distances of the restraints",
          command::Parameter( "filename", "\tfilename containing a list of pdb filenames", "pdb.ls")
        )
      );
      // end
      return s_flag;
    }

    //! @brief flag for specifying the number of desired restraints as fraction of number of residues in sses
    //! @return flag for specifying the number of desired restraints as fraction of number of residues in sses
    const util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagDataSetSizeFractionOfPool()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "data_set_size_fraction_of_sse_resis",
          "\tThe number of restraint should be based on the given fraction. The fraction is the fraction of restraints"
          " per residue in sses, where the residues in sses are given by the pool. If the pool contains overlapping "
          "sse definitions, a random non-overlapping subset is used for the number of restraints calculation. Don't"
          " forget to provide a pool when using this flag or the values from the other size flag will be used.",
          command::Parameter
          (
            "fraction", "\tThe fraction of restraints per residue in sses", "0.2"
          )
        )
      );

      // end
      return s_flag;
    }

    //! @brief flag for specifying the number of desired restraints as connecting all sses a given number of times
    //! @return flag for specifying the number of desired restraints as connecting all sses a given number of times
    const util::ShPtr< command::FlagInterface> &OptimizeDataSetPairwise::GetFlagDataSetSizeConnectSSEsInPool()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "data_set_size_connect_sses",
          "\tThe number of restraints should be enough to connect sses in the pool with each other the given number of"
          " times. If the pool contains overlapping "
          "sse definitions, a random non-overlapping subset is used for the number of restraints calculation. Don't"
          " forget to provide a pool when using this flag or the values from the other size flag will be used.",
          command::Parameter
          (
            "connection times", "\tnumber of times each sse should be connected with every other", "1"
          )
        )
      );

      // end
      return s_flag;
    }

    const ApplicationType OptimizeDataSetPairwise::OptimizeDataSetPairwise_Instance
    (
      GetAppGroups().AddAppToGroup( new OptimizeDataSetPairwise(), GetAppGroups().e_Restraint)
    );
  } // namespace app
} // namespace bcl
