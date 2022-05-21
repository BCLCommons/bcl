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

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_flex_field.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "chemistry/bcl_chemistry_fragment_evolve_implementations.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_mutate_interface.h"
#include "chemistry/bcl_chemistry_fragment_react.h"
#include "chemistry/bcl_chemistry_fragment_split_gadd_fragments.h"
#include "chemistry/bcl_chemistry_reaction_search.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_template_instantiations.h"
#include "pdb/bcl_pdb_factory.h"
#include "random/bcl_random_uniform_distribution.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_stopwatch.h"

// External includes - sorted alphabetically
#include "descriptor/bcl_descriptor_molecule_druglike.h"
#include <cstdio>
#include <iomanip>
namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MolInfo
    //! @brief A convenience class for holding molecule information
    //!
    //! @author geanesar, brownbp1
    //! @date 05/12/2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct MolInfo
    {

      // data
      std::string m_Identifier;
      chemistry::FragmentComplete m_Molecule;
      float m_Fitness;
      std::string m_History;
      size_t m_Age;

      // brief constructor
      MolInfo() :
        m_Identifier(),
        m_Molecule(),
        m_Fitness(),
        m_History(),
        m_Age( 0)
      {
      }

      // brief constructor with arguments
      MolInfo
      (
        const std::string &IDENTIFIER,
        const chemistry::FragmentComplete &MOLECULE,
        const float &FITNESS,
        const std::string &HISTORY = "",
        const size_t &AGE = 0
      ) :
        m_Identifier( IDENTIFIER),
        m_Molecule( MOLECULE),
        m_Fitness( FITNESS),
        m_History( HISTORY),
        m_Age( AGE)
      {
      }

      //! @brief less-than operator for MolInfos
      //! @return true if fitness of left operand is less than fitness of right operand
      bool operator <( const MolInfo &SECOND) const
      {
        return m_Fitness < SECOND.m_Fitness;
      }

      //! @brief greater-than operator for MolInfos
      //! @return true if fitness of left operand is greater than fitness of right operand
      bool operator >( const MolInfo &SECOND) const
      {
        return m_Fitness > SECOND.m_Fitness;
      }
    };

    // Forward declarations of comparison functions

    bool operator <( const util::ShPtr< MolInfo> &FIRST, const util::ShPtr< MolInfo> &SECOND);
    bool operator >( const util::ShPtr< MolInfo> &FIRST, const util::ShPtr< MolInfo> &SECOND);
    bool operator <( const std::pair< const MolInfo *, size_t> &FIRST, const std::pair< const MolInfo *, size_t> &SECOND);
    bool operator >( const std::pair< const MolInfo *, size_t> &FIRST, const std::pair< const MolInfo *, size_t> &SECOND);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EvoGen
    //! @brief A de-novo design application for small molecules using a stochastic search algorithm
    //!
    //! @author geanesar
    //! @date 05/12/2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class EvoGen :
      public Interface
    {

    protected:

      //! native target ligands
      mutable chemistry::FragmentEnsemble m_NativeLigands;

      //! native target pockets
      mutable chemistry::FragmentEnsemble m_BindingPocket;

    private:

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class MolOptimizer
      //! @brief The class that actually does the heavy-lifing part of the stochastic molecular structure search
      //!
      //! @author geanesar
      //! @date 11/11/2016
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class MolOptimizer
      {
      private:

        //! @brief methods for how molecules should be selected
        enum SelectionType
        {
          e_Top,                  //!< Select the top members
          e_Tournament,           //!< Select members by tournament
          s_NumberSelectionTypes
        };

        //! @brief how models are read in (BCL internal or external script)
        enum ModelType
        {
          e_Internal,         //!< BCL models only
          e_External,         //!< Launch an external script
          s_NumberModelTypes
        };

        //! @brief how parent molecules are treated
        enum ParentRetirementType
        {
          e_DoNotRetire,              //!< Retirement policy "none", keep all parents
          e_RetireProbabilisticByAge, //!< Retirement policy "probabilistic", retire parents by age
          e_RetireAll,                //!< Retirement policy "all", do not keep parents
          s_NumberRetirementTypes
        };

        //! @brief how the molecule evolution is balanced
        //! InsertMols always steady at 15%
        enum EvolutionBalance
        {
          e_AlchemicalMutate,
          e_ReactionInsertionOnly,
          e_ReactionDominant,               //!< 75% Reaction, 10% Recombination
          e_RecombinationDominant,          //!< 10% Reaction, 75% Recombination
          e_RecombinationInsertionOnly,
          e_ReactionRecombinationBalanced,  //!< 42.5% Reaction, 42.5% Recombination
          s_NumberEvolutionBalances
        };

        //! array for storing population information
        //! should contain m_FinalPopSize MolInfos per element
        std::vector< std::vector< MolInfo> > m_Populations;

        //! the maximum number of molecules to generate, per population
        size_t m_GenerateMax;

        //! the number of molecules to keep, per population
        size_t m_FinalPopSize;

        //! The reaction operation class, used for modifying structures
        chemistry::FragmentReact m_ReactOp;

        //! The alchemical mutate to use to combine molecules
        util::Implementation< chemistry::FragmentMutateInterface> m_Mutate; //!< obtains a implementation

        //! vector of molecules for addition/insertion operations
        util::ShPtrVector< chemistry::FragmentComplete> m_InsertMols;

        //! molecule scorer, if model type uses BCL-only models
        descriptor::CheminfoProperty m_Scorer;

        //! coordinates for the center of the binding pocket
        mutable linal::Vector3D m_ReferenceCoordinates;

        //! tracker for how many times different operations are used
        mutable std::vector< size_t> m_Operations;

        //! how molecules should be selected
        //! @see SelectionType
        SelectionType m_SelectionType;

        //! size of tournaments for the replacement step, if using tournament selection
        float m_ReplacementTournSizeFactor;

        //! size of the tournaments for each modification step, if using tournament selection
        float m_ModifyTournSizeFactor;

        //! the type of models that will be used
        //! @see ModelType
        ModelType m_ModelType;

        //! the shell command that should be run if using external scoring functions
        std::string m_ModelCmd;

        //! parent retirement policy
        //! @see ParentRetirementType
        ParentRetirementType m_RetirementType;

        //! how molecule evolution processes are balanced
        //! @see EvolutionBalance
        EvolutionBalance m_EvolutionBalanceType;

        //! vector of molecules for recombination operations
        mutable chemistry::FragmentEvolveImplementations m_RecombineOp;

        //! filename for where EvoGen activity should be logged
        std::string m_LogFile;

        //! output stream for EvoGen logging information
        io::OFStream m_LogStream;

      public:

        //! @brief constructor
        MolOptimizer() :
          m_GenerateMax( 10),
          m_FinalPopSize( 10),
          m_Operations( 3, size_t( 0)),
          m_SelectionType( e_Top),
          m_ModelType( e_Internal),
          m_ModelCmd(),
          m_RetirementType( e_DoNotRetire),
          m_EvolutionBalanceType( e_ReactionDominant),
          m_RecombineOp( chemistry::FragmentEvolveImplementations::EvolveType::e_Combine, false)
        {
        }

        //! @brief set up the initial population from an ensemble of molecules
        //! @details molecules must have >= 1 atoms or they will be discarded.
        //! @param MOLS the structures to use, will be copied
        void SetInitialPopulation( const chemistry::FragmentEnsemble &MOLS)
        {
          // check necessary components are available
          BCL_Assert( m_Scorer.IsDefined() || !m_ModelCmd.empty(), "Must set score descriptor before calling SetInitialPopulation");

          // clear all data and start with a fresh, empty, population
          m_Populations.clear();
          m_Populations.push_back( std::vector< MolInfo>());
          std::vector< MolInfo> &new_pop( m_Populations[ 0]);

          // perform checks on each molecule and add to new population with proper information
          size_t mol_no( 0);
          for
          (
            chemistry::FragmentEnsemble::const_iterator itr_mol( MOLS.Begin()), itr_mol_end( MOLS.End());
            itr_mol != itr_mol_end;
            ++itr_mol
          )
          {
            if( itr_mol->GetNumberAtoms())
            {
              new_pop.push_back
              (
                MolInfo( "P" + util::Format()( 0) + "M" + util::Format()( mol_no++), *itr_mol, 0, "Begin")
              );
            }
          }

          ScoreMolecules( new_pop);

          // sort molecules by fitness low->high
          std::sort( new_pop.begin(), new_pop.end());
        }

        //! @brief set filename for EvoGen log file
        //! @details this file is a json-formatted log file containing information about the generated molecules
        void SetLogFile( const std::string &FILENAME)
        {
          m_LogFile = FILENAME;
        }

        //! @brief open log file for writing; continues if file cannot be opened
        void StartLogging()
        {
          if( !io::File::TryOpenOFStream( m_LogStream, m_LogFile))
          {
            BCL_MessageStd( "Could not open " + m_LogFile + " for writing");
          }
        }

        //! @brief close/flush logging file stream
        void StopLogging()
        {
          if( m_LogStream.is_open())
          {
            io::File::CloseClearFStream( m_LogStream);
          }
        }

        //! @brief set molecule selection type to keep highest-scoring molecules
        void SetSelectionTop()
        {
          m_SelectionType = e_Top;
          BCL_MessageStd( "Selection type is: " + util::Format()( m_SelectionType));
        }

        //! @brief Set molecule replacement method to use tournament selection
        //! @param FACTOR the percentage of the available data that will be used in a single tournament round
        //! @details will warn if FACTOR < 0 or FACTOR > 1 and set FACTOR to 0.0 or 1.0, respectively.  The number of
        //!  molecules selected will be 1 <= NUMBER <= data size regardless of FACTOR
        void SetReplacementTypeTournament( const float &FACTOR)
        {
          float factor( FACTOR);
          if( FACTOR < 0.0)
          {
            BCL_MessageStd( "Tournament size factor is below 0.0, setting to 0");
            factor = 0.0;
          }
          else if( FACTOR > 1.0)
          {
            BCL_MessageStd( "Tournament size factor is above 1.0, setting to 1");
            factor = 1.0;
          }

          m_ReplacementTournSizeFactor = factor;
          m_SelectionType = e_Tournament;
          BCL_MessageStd( "Selection type is: " + util::Format()( m_SelectionType));
        }

        //! @brief Set molecule modification method to use tournament selection
        //! @param FACTOR the percentage of the available data that will be used in a single tournament round
        //! @details will warn if FACTOR < 0 or FACTOR > 1 and set FACTOR to 0.0 or 1.0, respectively.  The number of
        //!  molecules selected will be 1 <= NUMBER <= data size regardless of FACTOR
        void SetModifyTypeTournament( const float &FACTOR)
        {
          float factor( FACTOR);
          if( FACTOR < 0.0)
          {
            BCL_MessageStd( "Tournament size factor is below 0.0, setting to 0");
            factor = 0.0;
          }
          else if( FACTOR > 1.0)
          {
            BCL_MessageStd( "Tournament size factor is above 1.0, setting to 1");
            factor = 1.0;
          }
          m_ModifyTournSizeFactor = factor;
          m_SelectionType = e_Tournament;
        }

        //! @brief set molecule evolution type to primarily perform reactions
        void SetEvolutionBalanceAlchemicalMutate()
        {
          m_EvolutionBalanceType = e_AlchemicalMutate;
          BCL_MessageStd( "Evolution balance type is: AlchemicalMutate");
        }

        void SetEvolutionBalanceReactionDominant()
        {
          m_EvolutionBalanceType = e_ReactionDominant;
          BCL_MessageStd( "Evolution balance type is: ReactionDominant");
        }

        //! @brief set molecule evolution type to primarily perform reactions and insertions
        void SetEvolutionBalanceReactionInsertionOnly()
        {
          m_EvolutionBalanceType = e_ReactionInsertionOnly;
          BCL_MessageStd( "Evolution balance type is: ReactionInsertionOnly");
        }

        //! @brief set molecule evolution type to primarily perform recombinations
        void SetEvolutionBalanceRecombinationDominant()
        {
          m_EvolutionBalanceType = e_RecombinationDominant;
          BCL_MessageStd( "Evolution balance type is: RecombinationDominant");
        }

        //! @brief set molecule evolution type to primarily perform recombinations and insertions
        void SetEvolutionBalanceRecombinationInsertionOnly()
        {
          m_EvolutionBalanceType = e_RecombinationInsertionOnly;
          BCL_MessageStd( "Evolution balance type is: RecombinationInsertionOnly");
        }

        //! @brief set molecule evolution type to primarily perform balanced processes
        void SetEvolutionBalancedBalanced()
        {
          m_EvolutionBalanceType = e_ReactionRecombinationBalanced;
          BCL_MessageStd( "Evolution balance type is: ReactionRecombinationBalanced");
        }

        //! @brief use BCL-internal models for molecule scoring
        void SetModelTypeInternal()
        {
          m_ModelType = e_Internal;
        }

        //! @brief use external script for molecule scoring
        void SetModelTypeExternal()
        {
          m_ModelType = e_External;
        }

        //! @brief get the internal molecule data
        //! @return a vector containing population data for each iteration so far
        const std::vector< std::vector< MolInfo> > &GetMolInfos() const
        {
          return m_Populations;
        }

        //! @brief set the final size of each population.
        //! @param POP_SIZE the final size of populations
        void SetFinalPopSize( const size_t &POP_SIZE)
        {
          m_FinalPopSize = POP_SIZE;
        }

        //! @brief set the maximum number of molecules to generate during each iteration.  This will be pruned to
        //!  m_FinalPopSize afterwards
        //! @param MAX maximum number of molecules to generate
        void SetMaxToGenerate( const size_t &MAX)
        {
          m_GenerateMax = MAX;
        }

        //! @brief initialize the reaction operation structure
        //! @brief REACTANT_FILENAME filename from which to read reactant molecules
        //! @brief REACTION_DIRNAME directory in which RXN files should be found
        void SetupReactOperation( const std::string &REACTANT_FILENAME, const std::string &REACTION_DIRNAME)
        {
          m_ReactOp = chemistry::FragmentReact( chemistry::ReactionSearch( REACTANT_FILENAME, REACTION_DIRNAME));
        }

        //! @brief initialize the one-shot reaction operation structure
        //! @brief IMPLEMENTATION alchemical mutate implementation
        void SetupAlchemicalMutate( const std::string &IMPLEMENTATION)
        {
          // load fragment pool
          m_Mutate = IMPLEMENTATION;
        }

        //! @brief set the descriptor to use for scoring, if internal models are used
        //! @param DESCRIPTOR a string which should property encode a descriptor::CheminfoProperty
        //! @details this asserts that internal models are used to prevent programming errors
        void SetModelDescriptor( const std::string &DESCRIPTOR)
        {
          BCL_Assert
          (
            m_ModelType == e_Internal,
            "Programming error, tried setting model descriptor when not using BCL-internal models!"
          );

          m_Scorer.AssertRead( DESCRIPTOR);
        }

        //! @brief set the external script path if external scoring is used
        //! @param CMD the command to use.  Should be executed with the format: CMD <sdf from BCL> <output SDF>
        //! @details this asserts that internal models are used to prevent programming errors
        void SetModelCmd( const std::string &CMD)
        {
          BCL_Assert
          (
            m_ModelType == e_External,
            "Programming error, tried setting a model script when using BCL-internal models!"
          );
          m_ModelCmd = CMD;
        }

        //! @brief set retirement policy so that no parents are discarded
        void SetRetirementTypeNone()
        {
          m_RetirementType = e_DoNotRetire;
        }

        //! @brief set retirement policy so that parents are discarded probabilistically based on age
        void SetRetirementTypeProbabilistic()
        {
          m_RetirementType = e_RetireProbabilisticByAge;
        }

        //! @brief set retirement policy so that all parents are discarded
        void SetRetirementTypeAll()
        {
          m_RetirementType = e_RetireAll;
        }

        //! @brief set up the molecule insertion object
        //! @param INSERT_MOL_FILENAME sdf files which should be added to populations during a run
        void SetupInsertOperation( const std::string &INSERT_MOL_FILENAME)
        {
          io::IFStream in;
          io::File::MustOpenIFStream( in, INSERT_MOL_FILENAME);
          chemistry::FragmentFeed ff( in, sdf::e_Saturate);
          m_InsertMols.Reset();
          for( ; ff.NotAtEnd(); ++ff)
          {
            m_InsertMols.PushBack( util::CloneToShPtr( *ff));
          }
          io::File::CloseClearFStream( in);
        }

        //! @brief select a molecule by tournament selection
        //! @param MOLS the population to select from
        //! @param TOURN_SIZE the tournament size (0 <= TOURN_SIZE <= MOLS.size())
        //! @param IGNORE_INDICES the indices to ignore when selecting
        //! @return an index in the MOLS vector, or MOLS.size() if an error occurs
        //! @details if TOURN_SIZE == 0 then the highest-scoring molecule will be used
        size_t SelectMoleculeTournament
        (
          const std::vector< MolInfo> &MOLS,
          size_t TOURN_SIZE,
          const std::set< int> &IGNORE_INDICES = std::set< int>()
        ) const
        {
          if( MOLS.empty() || IGNORE_INDICES.size() == MOLS.size())
          {
            return MOLS.size();
          }

          // tournament size of 0 means no tournament, select best mol
          TOURN_SIZE = TOURN_SIZE > 0 ? std::min< size_t>( MOLS.size(), TOURN_SIZE) : MOLS.size();

          // vector of available indices
          std::vector< size_t> inds;
          inds.reserve( MOLS.size());

          // add all indices except those in IGNORE_INDICES
          for( size_t i( 0); i < MOLS.size(); ++i)
          {
            if( IGNORE_INDICES.find( i) == IGNORE_INDICES.end())
            {
              inds.push_back( i);
            }
          }

          // shuffle indices randomly and take a subset (if sufficient indices are around)
          if( TOURN_SIZE < MOLS.size())
          {
            std::random_shuffle( inds.begin(), inds.end());
            inds.resize( TOURN_SIZE);
          }

          // make a vector of MolInfo* and save the associated index in MOLS so that we can sort them by fitness
          // track fitness sums and minimum to properly calculate probabilities
          std::vector< std::pair< const MolInfo *, size_t> > m_ptrs;
          double total_fitness( 0);
          for( size_t i( 0); i < inds.size(); ++i)
          {
            m_ptrs.push_back( std::pair< const MolInfo *, size_t>( &MOLS[ inds[ i]], inds[ i]));
            total_fitness += MOLS[ inds[ i]].m_Fitness;
          }

          // sort the pointers from high to low by fitness
          std::sort( m_ptrs.begin(), m_ptrs.end(), std::greater< std::pair< const MolInfo *, size_t> >());

          // make sure total fitness is not too small for division purposes
          total_fitness = std::max< double>( 0.001, total_fitness);

          // probability of selecting the first member is proportional to its fitness relative to the whole data set,
          // but should not be lower than 10%
          double p( std::max< double>( 0.1, m_ptrs[ 0].first->m_Fitness / total_fitness));
          double one_minus_p( 1 - p);

          // test determines which element is selected, accum decides when test has picked the right element
          // the first element will have a probability of getting picked proportional to its fitness relative to all others
          // but should have at least a 10 percent chance of being selected

          double test( random::GetGlobalRandom().Random< double>( 0.0, 1.0)); // threshold value
          double accum( p); // accumulation of probability

          // the first element is selected with probability P, second with P*(1-P), third has P*(1-P)*(1-P), ...
          size_t n( 0), last_ptr( m_ptrs.size() - 1);
          for( ; test > accum && n < last_ptr; ++n, accum += accum * one_minus_p);

          return m_ptrs[ n].second;
        }

        //! @brief remove whitespace (via isspace) from a string
        //! @param STR the string to remove whitespace from
        //! @return STR without any whitespace
        std::string RemoveWhitespace( const std::string &STR) const
        {
          std::string str;
          str.reserve( STR.length());
          for
          (
            std::string::const_iterator itr( STR.begin()), itr_end( STR.end());
            itr != itr_end;
            ++itr
          )
          {
            if( !isspace( *itr))
            {
              str.push_back( *itr);
            }
          }
          return str;
        }

        //! @brief prepare a string for writing to CSV by escaping quotes
        std::string PrepareForCSV( const std::string &STR) const
        {
          std::string str( STR);
          size_t p = 0;
          bool needs_quotes( false);

          while( ( p = str.find( "\"", p)) != std::string::npos)
          {
            needs_quotes = true;
            str.insert( p, "\\");
          }

          // add quotes around
          if( str.find( ",") != std::string::npos)
          {
            needs_quotes = true;
          }

          if( needs_quotes)
          {
            str.insert( 0, "\"");
            str.append( "\"");
          }
          return str;
        }

        //! @brief escape quotes with '\' in a string
        //! @param STR the string to escape
        //! @return a copy of STR with escaped quotes
        std::string EscapeQuotes( const std::string &STR) const
        {
          std::string str( STR);
          size_t p = 0;
          while( ( p = str.find( "\"")) != std::string::npos)
          {
            str.insert( p, "\\");
          }
          while( ( p = str.find( "'")) != std::string::npos)
          {
            str.insert( p, "\\");
          }
          return str;
        }

        //! @brief execute a reaction/addition operation to generate new molecules
        //! @brief MOLS the vector to add the new molecules to
        //! @details this only generates structures and sets history values appropriately, but does nothing for ensuring
        //!  molecules are suitable for scoring or that their
        void GenerateMolecules( std::vector< MolInfo> &MOLS) const
        {
          if( m_Populations.empty())
          {
            return;
          }

          // reference to most recent population, should be sorted by score l->h
          const std::vector< MolInfo> &last_pop( m_Populations.back());
          if( last_pop.empty())
          {
            return;
          }

          //Select evolutionary process probabilities
          bool alchemy( false);
          float react( 0.0), recombine( 0.0);
          switch( m_EvolutionBalanceType)
          {
            case e_AlchemicalMutate:
            alchemy = true;
            break;
            case e_ReactionDominant:
            react = 0.75;
            recombine = 0.85;
            break;
            case e_ReactionInsertionOnly:
            react = 0.75;
            recombine = 0.75 + math::GetLowestBoundedValue< float>();
            break;
            case e_RecombinationDominant:
            react = 0.10;
            recombine = 0.85;
            break;
            case e_RecombinationInsertionOnly:
            react = math::GetLowestBoundedValue< float>();
            recombine = 0.85;
            break;
            case e_ReactionRecombinationBalanced:
            react = 0.425;
            recombine = 0.85;
            break;
            default:
            BCL_Exit( "No proper evolution balancing type was specified, exiting...", -1);
          }

          float op_test( random::GetGlobalRandom().Random< float>( 0, 1.0));
          size_t tourn_size( std::max< size_t>( 1, std::min< size_t>( m_ModifyTournSizeFactor * last_pop.size(), last_pop.size())));

          // one-shot reaction couplings using undefined atoms between two reagents via AddMedChem
          if( alchemy)
          {
            // probably needs an operations count

            // select molecule via tournament selection
            size_t mol_no( SelectMoleculeTournament( last_pop, tourn_size));

            // react the molecule with something random
            const chemistry::FragmentComplete &picked_mol( last_pop[ mol_no].m_Molecule);
//            m_Mutate->SetMutableElementsAtomIndices
//            (
//              picked_mol,
//              storage::Vector< chemistry::ElementType>(1, chemistry::GetElementTypes().e_Undefined)
//            );

            // mutate
            auto product( ( *m_Mutate)( picked_mol));

            // update info
            if( product.GetArgument().IsDefined())
            {
              std::stringstream hist;
              hist << "AlchemicalMutate,";
              hist << last_pop[ mol_no].m_Identifier;
              MOLS.push_back( MolInfo( "", *( product.GetArgument()), 0, hist.str()));
            }
          }

          // rxn-type reaction
          else if( op_test < react)
          {
            ++m_Operations[ 0];

            // select molecule via tournament selection
            size_t mol_no( SelectMoleculeTournament( last_pop, tourn_size));

            // react the molecule with something random
            const chemistry::FragmentComplete &picked_mol( last_pop[ mol_no].m_Molecule);
            storage::Pair< util::SiPtr< const chemistry::ReactionComplete>, chemistry::FragmentEnsemble> res
            (
              m_ReactOp.ReactRandom( picked_mol)
            );

            // a reaction may have generated more than one molecule, add each of them to the new population
            for
            (
              chemistry::FragmentEnsemble::const_iterator itr_mol( res.Second().Begin()), itr_mol_end( res.Second().End());
                itr_mol != itr_mol_end;
                ++itr_mol
            )
            {
              std::stringstream hist;
              hist << "React,";
              hist << RemoveWhitespace( res.First()->GetDescription());
              hist << "," << last_pop[ mol_no].m_Identifier;
              MOLS.push_back( MolInfo( "", *itr_mol, 0, hist.str()));
            }
          }
          // recombination
          else if( op_test < recombine)
          {
            ++m_Operations[ 1];

            // select molecule from tournament selection
            size_t mol_one( SelectMoleculeTournament( last_pop, tourn_size));
            const chemistry::FragmentComplete &picked_mol_one( last_pop[ mol_one].m_Molecule);
            util::ShPtrVector< chemistry::FragmentComplete> new_mols;

            // select recombination type randomly
            size_t recomb_type( random::GetGlobalRandom().Random< size_t>( 0, 2));
            if( recomb_type == 0 || recomb_type == 1)
            {
              // combine selected molecule with random molecule from the inserted (or "migrant") mols population
              m_RecombineOp = chemistry::FragmentEvolveImplementations( chemistry::FragmentEvolveImplementations::e_Combine, false);
              size_t mol_two( random::GetGlobalRandom().Random< size_t>( 0, m_InsertMols.GetSize() - 1));
              const chemistry::FragmentComplete &picked_mol_two( *m_InsertMols( mol_two));
              new_mols = m_RecombineOp.Combine( picked_mol_one, picked_mol_two, 50);
            }
            else if( recomb_type == 3)
            {
              // split GADD fragments from a random inserted (or "migrant) molecule and
              m_RecombineOp = chemistry::FragmentEvolveImplementations( chemistry::FragmentEvolveImplementations::e_FragAdd, false);
//              size_t mol_two( random::GetGlobalRandom().Random< size_t>( 0, m_InsertMols.GetSize() - 1));
//              const chemistry::FragmentComplete &picked_mol_two( *m_InsertMols( mol_two));

              //List of component fragment indices
              storage::List< storage::Vector< size_t> > gadd_fragment_indices;

              // GADD fragment rings
              if( random::GetGlobalRandom().Random< size_t>( 0, 1) == 0)
              {
                chemistry::FragmentSplitGADDFragments splitter;
//                splitter.GetComponentVertices( picked_mol_two, picked_mol_two.)
              }

              // GADD fragment chains
              else
              {
                chemistry::FragmentSplitGADDFragments splitter;
              }

//              new_mols = m_RecombineOp.MutateAdd()

            }
            else if( recomb_type == 2)
            {
              m_RecombineOp = chemistry::FragmentEvolveImplementations( chemistry::FragmentEvolveImplementations::e_FragDel, false);
              new_mols = m_RecombineOp.MutateDel( picked_mol_one);
            }

            // now add the new compounds to the mix
            for
            (
              util::ShPtrVector< chemistry::FragmentComplete>::iterator new_mols_itr( new_mols.Begin()), new_mols_itr_end( new_mols.End());
                new_mols_itr != new_mols_itr_end;
                ++new_mols_itr
            )
            {
              std::stringstream hist;
              hist << "Recombine,";
              hist << last_pop[ mol_one].m_Identifier;
//              hist << "," << m_InsertMols(mol_two)->GetName();
              MOLS.push_back( MolInfo( "", **new_mols_itr, 0, hist.str()));
            }
          }

          else
          {
            ++m_Operations[ 2];
            if( m_InsertMols.GetSize() > 1)
            {
              size_t mol_no( random::GetGlobalRandom().Random< size_t>( 0, m_InsertMols.GetSize() - 1));
              MOLS.push_back( MolInfo());
              MolInfo &mol = MOLS.back();
              mol.m_Molecule = *m_InsertMols( mol_no);
              mol.m_History = "Insert,M" + util::Format()( mol_no);
            }
            else if( m_InsertMols.GetSize() == 1)
            {
              MOLS.push_back( MolInfo());
              MolInfo &mol = MOLS.back();
              mol.m_Molecule = *m_InsertMols( 0);
              mol.m_History = "Insert,M" + util::Format()( 0);
            } // else noop
          }
        }

        //! @brief helper function for scoring a single molecule using internal BCL models
        //! @param MOL the molecule to score
        //! @return the score of the molecule
        //! @details does not do any preprocessing of the molecule, this is up to the caller
        float ScoreMoleculeInternal( const chemistry::FragmentComplete &MOL) const
        {
          BCL_Assert( m_Scorer.IsDefined(), "scorer was not defined");
          return m_Scorer->SumOverObject( MOL)( 0);
        }

        //! @brief helper function for scoring a single molecule using an external script
        //! @param MOL_INFOS the molecules to score
        //! @details this will update the MOL_INFOS structure directly with the new fitness scores as read from the
        //!  output of the called script
        void ScoreMoleculesExternal( std::vector< MolInfo> &MOL_INFOS) const
        {
          BCL_Assert( !m_ModelCmd.empty(), "No model command given!");

          time_t cur_time;
          time( &cur_time);
          size_t rand_no( random::GetGlobalRandom().Random< size_t>( 0, 2 << 28));

          // try to make a sufficiently random filename
          std::string file_basename( "evogen_score_" + util::Format()( cur_time) + util::Format()( rand_no));

          std::string filename_out = "/tmp/" + file_basename + "_out.sdf";
          std::string filename_in = "/tmp/" + file_basename + "_in.sdf";

          // keep track of molecule names in case they get rearranged
          std::map< std::string, MolInfo *> mol_names;

          // write out molecule data to sdf
          io::OFStream out;
          io::File::MustOpenOFStream( out, filename_out);
          size_t mol_no( 0);
          for
          (
            std::vector< MolInfo>::iterator itr_mi( MOL_INFOS.begin()), itr_mi_end( MOL_INFOS.end());
            itr_mi != itr_mi_end;
            ++itr_mi, ++mol_no
          )
          {
            BCL_Assert
            (
              !itr_mi->m_Identifier.empty(),
              "Programming error: molecule " + util::Format()( mol_no) + " has an empty identifier!"
            );

            // remove stale properties
            if( !itr_mi->m_Molecule.IsPropertyStored( "TmpEvoGenIdentifier"))
            {
              itr_mi->m_Molecule.RemoveProperty( "TmpEvoGenIdentifier");
            }

            // add an identifier to the molecule
            itr_mi->m_Molecule.StoreProperty( "TmpEvoGenIdentifier", itr_mi->m_Identifier);

            // All identifiers should be unique so this should never actually happen, but double check
            if( mol_names.find( itr_mi->m_Identifier) != mol_names.end())
            {
              BCL_MessageStd( "Multiple molecules were given with ID \"" + itr_mi->m_Identifier + "\"!");
            }
            mol_names[ itr_mi->m_Identifier] = &( *itr_mi);

            itr_mi->m_Molecule.WriteMDL( out);
          }

          io::File::CloseClearFStream( out);

          // execute the external command: cmd <infile to prog> <outfile from prog>
          std::string cmd( m_ModelCmd + " " + filename_out + " " + filename_in);
          int err = system( cmd.c_str());
          BCL_Assert( err == 0, "Could not execute command \"" + cmd + "\", exitted with status " + util::Format()( err));

          // read in what cmd wrote out
          io::IFStream in;
          io::File::MustOpenIFStream( in, filename_in);
          chemistry::FragmentEnsemble in_ens( in, sdf::HydrogenHandlingPref::e_Saturate);

          // check the input
          size_t ens_size( in_ens.GetSize());
          BCL_Assert
          (
            ens_size == MOL_INFOS.size(),
            "error reading output from command: output contained " + util::Format()( MOL_INFOS.size()) + " molecules "
            " but after scoring the input contained " + util::Format()( ens_size) + " molecules"
          );

          // determine the score for each molecule by checking for their identifiers
          mol_no = 0;
          std::map< std::string, MolInfo *>::iterator map_end( mol_names.end());
          for
          (
            chemistry::FragmentEnsemble::const_iterator itr_mol( in_ens.Begin()), itr_mol_end( in_ens.End());
            itr_mol != itr_mol_end;
            ++itr_mol, ++mol_no
          )
          {
            BCL_Assert
            (
              itr_mol->IsPropertyStored( "EvoGenFitness"),
              "Molecule " + util::Format()( mol_no) + " did not have an EvoGenFitness property"
            );
            BCL_Assert
            (
              itr_mol->IsPropertyStored( "TmpEvoGenIdentifier"),
              "Could not determine the identity of molecule " + util::Format()( mol_no) + " (no TmpEvoGenIdentifier property)"
            );

            std::string score_str( RemoveWhitespace( itr_mol->GetMDLProperty( "EvoGenFitness")));
            float score;

            BCL_Assert
            (
              util::TryConvertFromString( score, score_str, util::GetLogger()),
              "Could not convert EvoGenFitness with value \"" + score_str + "\" to float (for molecule " + util::Format()( mol_no) + ")"
            );

            std::string m_name( RemoveWhitespace( itr_mol->GetMDLProperty( "TmpEvoGenIdentifier")));
            std::map< std::string, MolInfo *>::iterator itr_map( mol_names.find( m_name));
            BCL_Assert
            (
              itr_map != map_end,
              "Could not find a molecule with identifier \"" + m_name + "\" in the MolInfo vector (out of " + util::Format()( mol_names.size()) + ")"
            );

            // update the molecule with whatever was provided by the scoring script
            itr_map->second->m_Molecule = *itr_mol;
            itr_map->second->m_Fitness = score;
            itr_map->second->m_Molecule.RemoveProperty( "TmpEvoGenIdentifier");
          }
          remove( filename_out.c_str());
          remove( filename_in.c_str());
        }

        //! @brief copy the highest scoring molecules from one population to a number
        //! @param FROM the population to copy from
        //! @param TO the population to copy to
        //! @param NUM the number to copy
        void CopyHighestScoring( std::vector< MolInfo> &FROM, std::vector< MolInfo> &TO, size_t NUM)
        {
          std::sort( FROM.begin(), FROM.end());
          size_t n_added( 0);
          for
          (
            std::vector< MolInfo>::const_reverse_iterator itr( FROM.rbegin()), itr_end( FROM.rend());
            itr != itr_end && n_added < NUM;
            ++itr, ++n_added
          )
          {
            TO.push_back( *itr);
          }
        }

        //! @brief picks molecules from FROM and copies them into TO
        //! @param FROM the population to copy molecules from
        //! @param TO the population to copy molecules to
        //! @param NUM the number to select
        //! @param TOURN_SIZE the tournament size to use
        //! @param REPLACEMENT whether to replace selected molecules (i.e. allow multiple picking)
        void CopyByTournament
        (
          const std::vector< MolInfo> &FROM,
          std::vector< MolInfo> &TO,
          size_t NUM,
          size_t TOURN_SIZE,
          const bool &REPLACEMENT = false
        ) const
        {
          size_t n_added( 0), n_left( FROM.size());
          std::set< int> already_picked;
          while( n_left > 0 && n_added < NUM)
          {
            size_t pick( SelectMoleculeTournament( FROM, TOURN_SIZE, already_picked));
            TO.push_back( FROM[ pick]);

            if( !REPLACEMENT)
            {
              already_picked.insert( pick);
              --n_left;
            }

            ++n_added;
          }
          std::sort( TO.begin(), TO.end());
        }

        //! @brief helper funciton to score a vector of MolInfos
        //! @param MOL_INFOS the MolInfos containing molecules to score
        //! @details this updates the MolInfo.m_Fitness function of all molecules
        void ScoreMolecules( std::vector< MolInfo> &MOL_INFOS)
        {
          switch( m_ModelType)
          {
            case e_Internal:
              for
              (
                std::vector< MolInfo>::iterator itr_mi( MOL_INFOS.begin()), itr_mi_end( MOL_INFOS.end());
                itr_mi != itr_mi_end;
                ++itr_mi
              )
              {
                itr_mi->m_Fitness = ScoreMoleculeInternal( itr_mi->m_Molecule);
              }
              break;
            case e_External:
              ScoreMoleculesExternal( MOL_INFOS);
              break;
            default:
              BCL_Exit( "Neither internal nor external model was specified", -1);
          }
        }

        //! @brief write data to the json log file
        //! @param STR the string to write
        void WriteLog( const std::string &STR)
        {
          if( m_LogStream.is_open())
          {
            m_LogStream << STR;
          }
        }

        //! @brief execute an iteration of molecule generation/scoring/replacement
        //! @return 0 on success, negative value on error
        int Next()
        {
          // druglikeness hack
          static descriptor::MoleculeDruglike druglike( false);

          std::vector< size_t> prev_ops( m_Operations);
          std::vector< MolInfo> next_pop;
          next_pop.reserve( m_GenerateMax + ( m_Populations.empty() ? 0 : m_Populations.back().size()));
          chemistry::ConstitutionSet unique_mols;

//          int tries = 0;
          int mol_no = 0;

          // do this until we reach the specified number of molecules
          while( next_pop.size() < m_GenerateMax)
          {

            std::vector< MolInfo> mols;
            GenerateMolecules( mols);

            // add the newly generated molecules to the population one by one...
            for( size_t i( 0); i < mols.size(); ++i)
            {

              util::ShPtr< chemistry::FragmentComplete> final_mol;
              m_Mutate.IsDefined() ?
              final_mol = util::ShPtr< chemistry::FragmentComplete> ( util::ShPtr< chemistry::FragmentComplete>( new chemistry::FragmentComplete( mols[i].m_Molecule))) :
              final_mol = util::ShPtr< chemistry::FragmentComplete> ( chemistry::FragmentEvolveBase::FinalizeMolecule( mols[ i].m_Molecule));
              if( final_mol.IsDefined())
              {
                mols[ i].m_Molecule = *final_mol;

                // make sure we don't generate duplicates in a single population
                if
                (
                    mols[ i].m_Molecule.GetNumberAtoms() &&
                    unique_mols.Insert( chemistry::FragmentConstitutionShared( mols[ i].m_Molecule)).second &&
                    druglike( mols[ i].m_Molecule)( 0)
                )
                {
                  mols[ i].m_Identifier = "P" + util::Format()( m_Populations.size()) + "M" + util::Format()( mol_no++);
                  next_pop.push_back( mols[ i]);
//                  tries = 0;
                }
              }
            }

            // if we are failing really hard, then we should probably just abort...
//            BCL_Assert( ++tries < 100, "Tried to generate a unique molecule, took more than 100 tries, bailing out");

            // print off a status message so the users don't get antsy
            if( next_pop.size() % 10)
            {
              util::GetLogger().LogStatus
              (
                "Generated " + util::Format()( next_pop.size() + 1) + " molecules for population "
                + util::Format()( m_Populations.size())
              );
            }
          }

          storage::Vector< storage::Triplet< chemistry::FragmentComplete, chemistry::FragmentComplete, double>> refined_mols( next_pop.size());

          // score all molecules as a batch
          BCL_MessageStd( "Scoring molecules");
          ScoreMolecules( next_pop);

          // Write information to the log file in json format.  This is ugly as sin and should be re-done in the future
          WriteLog( "{ \"pop_number\": " + util::Format()( m_Populations.size()) + ",\n");
          WriteLog( "  \"pre_replacement\": [\n");
          for( size_t i( 0); i < next_pop.size(); ++i)
          {
            WriteLog( "      { \"id\": \"" + next_pop[ i].m_Identifier + "\",\n");
            WriteLog( "        \"history\": \"" + EscapeQuotes( next_pop[ i].m_History) + "\",\n");
            WriteLog( "        \"age\": " + util::Format()( next_pop[ i].m_Age) + ",\n");
            WriteLog( "        \"score\": " + util::Format()( next_pop[ i].m_Fitness) + "}");
            if( i < ( next_pop.size() - 1))
            {
              WriteLog( ",\n");
            }
          }
          WriteLog( "  ],\n");

          // parent replacement: add members from the old population to the new population
          if( m_RetirementType != e_RetireAll)
          {
            BCL_MessageStd( "Copying parents to the new population");
            WriteLog( "  \"copied_parents\": [\n");

            const std::vector< MolInfo> &parent_pop( m_Populations.back());
            for( size_t i( 0); i < parent_pop.size(); ++i)
            {
              double r( 0.5); // test value for keeping a molecule; initial value is a dummy value for Retire None policy
              double cutoff( 1.0); // test threshold; r must be below this value to keep the molecule

              // if retirement policy is probabilistic then "roll the dice" to see if this molecule survives
              if( m_RetirementType == e_RetireProbabilisticByAge)
              {
                r = random::GetGlobalRandom().Random< double>( 0, 1.0); // assign a real test value
                cutoff = exp( double( -0.5) * parent_pop[ i].m_Age); // assign a real test threshold value based on age
              }

              // make sure the molecule passed the test, and it's a unique structure
              if( r < cutoff && unique_mols.Insert( chemistry::FragmentConstitutionShared( parent_pop[ i].m_Molecule)).second)
              {
                //BCL_MessageStd( "  Keeping parent " + parent_pop[ i].m_Identifier);
                next_pop.push_back( parent_pop[ i]);
                next_pop.back().m_Age++;

                // write to the log file
                WriteLog( "      { \"id\": \"" + parent_pop[ i].m_Identifier + "\",\n");
                WriteLog( "        \"history\": \"" + EscapeQuotes( parent_pop[ i].m_History) + "\",\n");
                WriteLog( "        \"age\": " + util::Format()( parent_pop[ i].m_Age) + ",\n");
                WriteLog( "        \"score\": " + util::Format()( parent_pop[ i].m_Fitness) + "}");
                if( i < ( parent_pop.size() - 1))
                {
                  WriteLog( ",\n");
                }

              }
              else
              {
                // status message to make users happy
                BCL_MessageStd( "  Retiring parent " + parent_pop[ i].m_Identifier);
              }
            }
            // additional closing braces to make sure json format is ok
            WriteLog( "\n  ],\n");
          }

          size_t new_pop_no( m_Populations.size()); // the new population number (same as current size since we add one)
          m_Populations.push_back( std::vector< MolInfo>()); // allocate space for new population
          m_Populations[ new_pop_no].reserve( m_FinalPopSize); // reserve space in the new pop for all mols

          BCL_MessageStd( "Downsampling the population");
          switch( m_SelectionType)
          {
            case e_Top:
              CopyHighestScoring( next_pop, m_Populations[ new_pop_no], m_FinalPopSize);
              //CopyByTournament( next_pop, m_Populations[ new_pop_no], m_FinalPopSize, next_pop.size());
              break;
            case e_Tournament:
              {
                size_t tourn_size( std::min< size_t>( std::max< size_t>( next_pop.size() * m_ReplacementTournSizeFactor, 1), next_pop.size()));
                CopyByTournament( next_pop, m_Populations[ new_pop_no], m_FinalPopSize, tourn_size);
              }
              break;
            default:
              return -1;
          }

          // warn the user if there aren't enough molecules
          if( m_Populations[ new_pop_no].size() < m_FinalPopSize)
          {
            BCL_MessageStd
            (
              "Final population size for " + util::Format()( new_pop_no) + " is not big enough (is " +
              util::Format()( m_Populations[ new_pop_no].size()) + ", should be " +
              util::Format()( m_FinalPopSize) + ")"
            );
          }
          else
          {
            BCL_MessageStd( "Population " + util::Format()( new_pop_no) + " contains " + util::Format()( m_Populations[ new_pop_no].size()) + " molecules");
          }

          // write more log information...
          WriteLog( "  \"post_replacement\": [\n");
          for( size_t i( 0); i < m_Populations[ new_pop_no].size(); ++i)
          {
            WriteLog( "      { \"id\": \"" + m_Populations[ new_pop_no][ i].m_Identifier + "\",\n");
            WriteLog( "        \"history\": \"" + EscapeQuotes( m_Populations[ new_pop_no][ i].m_History) + "\",\n");
            WriteLog( "        \"age\": " + util::Format()( m_Populations[ new_pop_no][ i].m_Age) + ",\n");
            WriteLog( "        \"score\": " + util::Format()( m_Populations[ new_pop_no][ i].m_Fitness) + "}");
            if( i < ( m_Populations[ new_pop_no].size() - 1))
            {
              WriteLog( ",\n");
            }
          }
          WriteLog( "\n  ]\n}\n");

          // tell the user how many reactions/insertions were done
          BCL_MessageStd( "Ran " + util::Format()( m_Operations[ 0] - prev_ops[ 0]) + " reactions");
          BCL_MessageStd( "Ran " + util::Format()( m_Operations[ 1] - prev_ops[ 1]) + " recombinations");
          BCL_MessageStd( "Ran " + util::Format()( m_Operations[ 2] - prev_ops[ 2]) + " inserts");
          return 0;
        }
      }; // end MolOptimizer

    //////////
    // data //
    //////////

      util::ShPtr< command::FlagInterface> m_StartingMolsFlag;
      util::ShPtr< command::FlagInterface> m_InsertMolsFlag;
      util::ShPtr< command::FlagInterface> m_ReactionsFlag;
      util::ShPtr< command::FlagInterface> m_ReagentsFlag;
      util::ShPtr< command::FlagInterface> m_ModelDirectoryFlag;
      util::ShPtr< command::FlagInterface> m_ModelPrefixFlag;
      util::ShPtr< command::FlagInterface> m_ExtScoringCommandFlag;
      util::ShPtr< command::FlagInterface> m_DruglikenessModelFlag;
      util::ShPtr< command::FlagInterface> m_WeightScoreFlag;
      util::ShPtr< command::FlagInterface> m_OutputDirFlag;
      util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;
      util::ShPtr< command::FlagInterface> m_PopulationSizeFlag;
      util::ShPtr< command::FlagInterface> m_OversampleFlag;
      util::ShPtr< command::FlagInterface> m_ShuffleFlag;
      util::ShPtr< command::FlagInterface> m_IterationsFlag;
      util::ShPtr< command::FlagInterface> m_ActiveCutoffFlag;
      util::ShPtr< command::FlagInterface> m_SelectionTypeFlag;
      util::ShPtr< command::FlagInterface> m_EvolutionBalanceTypeFlag;
      util::ShPtr< command::FlagInterface> m_ReplaceTournSizeFlag;
      util::ShPtr< command::FlagInterface> m_ModifyTournSizeFlag;
      util::ShPtr< command::FlagInterface> m_RetirementTypeFlag;
      util::ShPtr< command::FlagInterface> m_ImplementationFlag;

    public:

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      EvoGen();

      //! @brief clone function
      //! @return a pointer to a copy of this class
      EvoGen *Clone() const
      {
        return new EvoGen( *this);
      }

      //! @brief Get the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      //! @return a ShPtr containing the commands
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_StartingMolsFlag);
        sp_cmd->AddFlag( m_InsertMolsFlag);
        sp_cmd->AddFlag( m_ReactionsFlag);
        sp_cmd->AddFlag( m_ReagentsFlag);
        sp_cmd->AddFlag( m_ModelDirectoryFlag);
        sp_cmd->AddFlag( m_ModelPrefixFlag);
        sp_cmd->AddFlag( m_ExtScoringCommandFlag);
        sp_cmd->AddFlag( m_DruglikenessModelFlag);
        sp_cmd->AddFlag( m_WeightScoreFlag);
        sp_cmd->AddFlag( m_OutputDirFlag);
        sp_cmd->AddFlag( m_OutputPrefixFlag);
        sp_cmd->AddFlag( m_PopulationSizeFlag);
        sp_cmd->AddFlag( m_OversampleFlag);
        sp_cmd->AddFlag( m_ShuffleFlag);
        sp_cmd->AddFlag( m_IterationsFlag);
        sp_cmd->AddFlag( m_ActiveCutoffFlag);
        sp_cmd->AddFlag( m_SelectionTypeFlag);
        sp_cmd->AddFlag( m_EvolutionBalanceTypeFlag);
        sp_cmd->AddFlag( m_ReplaceTournSizeFlag);
        sp_cmd->AddFlag( m_ModifyTournSizeFlag);
        sp_cmd->AddFlag( m_RetirementTypeFlag);
        sp_cmd->AddFlag( m_ImplementationFlag);

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      std::string SetupDescriptorString() const
      {

        std::string descriptor_str( "Combine(Multiply(");

        std::string model_prefix( m_ModelPrefixFlag->GetFlag() ? m_ModelPrefixFlag->GetFirstParameter()->GetValue() : "model");
        descriptor_str += "PredictionMean(storage=File(directory="
                          + m_ModelDirectoryFlag->GetFirstParameter()->GetValue()
                          + ",prefix=" + model_prefix
                          +"))";
        if( m_WeightScoreFlag->GetFlag())
        {
          descriptor_str += ",WithinRangeSmooth(descriptor=Weight,begin=150,end=550,left width=50,right width=50)";
        }
        descriptor_str += "))";
        return descriptor_str;
      }

      void SetupInitialPopulation( chemistry::FragmentEnsemble &INIT_POP, const storage::Vector< chemistry::FragmentComplete> &MOLECULES, size_t POP_SIZE) const
      {
        for( size_t mol_num( 0), added( 0); mol_num < MOLECULES.GetSize() && added < POP_SIZE; ++mol_num)
        {
          if( MOLECULES( mol_num).GetNumberAtoms() > 0)
          {
            INIT_POP.PushBack( MOLECULES( mol_num));
            ++added;
          }
        }
      }

      //! @brief the main part of the application
      //! @return 0 on success
      int Main() const
      {

      /////////////////////////
      // Initial population  //
      /////////////////////////

        BCL_MessageStd( "Reading in molecules for starting population");

        //! Molecules that will be put into initial population
        storage::Vector< chemistry::FragmentComplete> molecules;

        //! @brief the population size as requested by user
        size_t requested_pop_size( m_PopulationSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>());

        // If shuffle is specified, read in everything.  Otherwise only read in however many the user wants in a population
        size_t n_to_read( m_ShuffleFlag->GetFlag() ? std::numeric_limits< size_t>::max() : requested_pop_size);

        // add indices so that the (alex stopped his comment here, indicating that there is no good reason to add indices - Ben)
        chemistry::FragmentFeed input_mols( m_StartingMolsFlag->GetStringList(), sdf::e_Saturate, n_to_read);
        for( size_t mol_no( 0); input_mols.NotAtEnd(); ++input_mols, ++mol_no)
        {
          if( input_mols->GetNumberAtoms() > 0)
          {
            chemistry::FragmentComplete next_mol( *input_mols);
            if( next_mol.IsPropertyStored( "EvoGenIdentifier"))
            {
              next_mol.RemoveProperty( "EvoGenIdentifier");
            }
            next_mol.StoreProperty( "EvoGenIdentifier", "P" + util::Format()( 0) + "M" + util::Format()( mol_no));
            molecules.PushBack( next_mol);
          }
        }

        BCL_MessageStd( "Read in a total of " + util::Format()( molecules.GetSize()) + " molecules");

        // If shuffling is specified, randomize the molecules
        if( m_ShuffleFlag->GetFlag())
        {
          BCL_MessageStd( "Shuffling molecules");
          molecules.Shuffle();
        }

        //! @brief the population size that will be used
        size_t set_pop_size = requested_pop_size;

        // Size of the initial population
        // Warn if there are not enough molecules to fill the initial population to user's request
        if( molecules.GetSize() < set_pop_size)
        {
          BCL_MessageCrt
          (
            "WARNING: population size specified as " + util::Format()( requested_pop_size)
            + " but input files only contain " + util::Format()( molecules.GetSize()) + "."
          );
        }

        chemistry::FragmentEnsemble init_pop;

        SetupInitialPopulation( init_pop, molecules, set_pop_size);

        BCL_MessageStd( "Initial population contains " + util::Format()( init_pop.GetSize()) + " members");

        // Write the status file to wherever the output filename will be written
        std::string output_dir( m_OutputDirFlag->GetFlag() ? m_OutputDirFlag->GetFirstParameter()->GetValue() : ".");

        std::string prefix( output_dir + "/");
        std::string out_prefix( m_OutputPrefixFlag->GetFirstParameter()->GetValue());
        std::string final_output_fn( prefix + out_prefix + "_final.sdf.gz");
        std::string checkpoint_filename( prefix + out_prefix + "_checkpoint.sdf.gz");

        std::string details_prefix( prefix + out_prefix + "_details");
        std::string member_history_filename = details_prefix + "_histories.txt";
        std::string member_details_filename = details_prefix + "_member_details.csv";
        std::string population_stats_filename = details_prefix + "_pop_stats.csv";
        std::string active_molecules_filename = details_prefix + "_active.sdf";

        std::string pop_sdf_prefix = details_prefix + "_sdf";

        util::ShPtr< chemistry::FragmentEnsemble> active_mols( new chemistry::FragmentEnsemble());

        BCL_MessageStd( "");
        BCL_MessageStd( "=============================");
        BCL_MessageStd( "RUN DETAILS:");
        BCL_MessageStd( "  Output file: --------------------- " + final_output_fn)
        BCL_MessageStd( "  Checkpoint file: ----------------- " + checkpoint_filename);

        if( !pop_sdf_prefix.empty())
        {
          BCL_MessageStd( "    (details file) Population SDFs:  " + pop_sdf_prefix + "_pop_X.sdf.gz");
        }

        BCL_MessageStd( "=============================\n");

        util::Stopwatch timer
        (
          "Molecule generation",
          util::Message::e_Standard,
          false
        );

        timer.Start();

        // Execute the evolutionary algorithm
        MolOptimizer optimizer;

        // set up the population filter/selector
        std::string sel_type( m_SelectionTypeFlag->GetFirstParameter()->GetValue());
        if( sel_type == "top")
        {
          optimizer.SetSelectionTop();
        }
        else if( sel_type == "tournament")
        {
          optimizer.SetReplacementTypeTournament( m_ReplaceTournSizeFlag->GetFirstParameter()->GetNumericalValue< float>());
          optimizer.SetModifyTypeTournament( m_ModifyTournSizeFlag->GetFirstParameter()->GetNumericalValue< float>());
        }
        else
        {
          BCL_MessageStd( "Unknown selection type \"" + m_SelectionTypeFlag->GetFirstParameter()->GetValue() + "\" specified");
          return -1;
        }

        // set up the evolution balance type
        std::string evbal_type( m_EvolutionBalanceTypeFlag->GetFirstParameter()->GetValue());
        if( evbal_type == "alchemical_mutate")
        {
          optimizer.SetEvolutionBalanceAlchemicalMutate();
        }
        else if( evbal_type == "reaction_dominant")
        {
          optimizer.SetEvolutionBalanceReactionDominant();
        }
        else if( evbal_type == "reaction_insertion_only")
        {
          optimizer.SetEvolutionBalanceReactionInsertionOnly();
        }
        else if( evbal_type == "recombination_dominant")
        {
          optimizer.SetEvolutionBalanceRecombinationDominant();
        }
        else if( evbal_type == "recombination_insertion_only")
        {
          optimizer.SetEvolutionBalanceRecombinationInsertionOnly();
        }
        else if( evbal_type == "balanced")
        {
          optimizer.SetEvolutionBalancedBalanced();
        }
        else
        {
            BCL_MessageStd( "Unknown evolution balance type \"" + m_EvolutionBalanceTypeFlag->GetFirstParameter()->GetValue() + "\" specified");
            return -1;
        }

        // determine retirement type
        std::string ret_type( m_RetirementTypeFlag->GetFirstParameter()->GetValue());
        if( ret_type == "all")
        {
          optimizer.SetRetirementTypeAll();
        }
        else if( ret_type == "probabilistic")
        {
          optimizer.SetRetirementTypeProbabilistic();
        }
        else if( ret_type == "none")
        {
          optimizer.SetRetirementTypeNone();
        }
        else
        {
          BCL_Exit( "Retirement type not recognized", -1);
        }

        // determine population sizes and sampling factors
        float oversample_factor( m_OversampleFlag->GetFirstParameter()->GetNumericalValue< float>());
        optimizer.SetFinalPopSize( set_pop_size);
        optimizer.SetMaxToGenerate( set_pop_size * oversample_factor);

        // set up scoring functions
        BCL_Assert
        (
          ( m_ExtScoringCommandFlag->GetFlag() || m_ModelDirectoryFlag->GetFlag())
            && !( m_ExtScoringCommandFlag->GetFlag() && m_ModelDirectoryFlag->GetFlag()),
          "Exactly one of " + m_ExtScoringCommandFlag->GetName() + " or " + m_ModelDirectoryFlag->GetName() + " must be specified"
        );

        // get the model directory
        if( m_ModelDirectoryFlag->GetFlag())
        {
          // set up fitness function
          std::string descriptor_str = SetupDescriptorString();
          BCL_MessageStd( "Score descriptor: " + descriptor_str);
          optimizer.SetModelTypeInternal();
          optimizer.SetModelDescriptor( descriptor_str);
        }
        else if( m_ExtScoringCommandFlag->GetFlag())
        {
          optimizer.SetModelTypeExternal();
          optimizer.SetModelCmd( m_ExtScoringCommandFlag->GetFirstParameter()->GetValue());
        }

        // set up reaction/insertion operations
        optimizer.SetupReactOperation( m_ReagentsFlag->GetFirstParameter()->GetValue(), m_ReactionsFlag->GetFirstParameter()->GetValue());
        optimizer.SetupInsertOperation( m_InsertMolsFlag->GetFirstParameter()->GetValue());
        optimizer.SetupAlchemicalMutate( m_ImplementationFlag->GetFirstParameter()->GetValue());

        // consider adding flag later to allow recombination
        // with a special subset of fragments

        // set the initial population
        optimizer.SetInitialPopulation( init_pop);

        // determine the number of iterations to run for
        size_t n_iters( m_IterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>());

        // write out initial population scores before we begin
        io::OFStream out;
        io::OFStream out_csv;
        io::File::MustOpenOFStream( out_csv, pop_sdf_prefix + "_all_scores.csv");

        out_csv << "pop_no,fitness\n";

        { // this is here on purpose so that we can re-use variable names later
          std::string out_filename( pop_sdf_prefix + "_0.sdf.gz");
          const std::vector< MolInfo> &mols( optimizer.GetMolInfos().back());

          io::File::MustOpenOFStream( out, out_filename);
          for( size_t i( 0), end_i( mols.size()); i < end_i; ++i)
          {
            chemistry::FragmentComplete mol( mols[ i].m_Molecule);
            if( mol.IsPropertyStored( "EvoGenFitness"))
            {
              mol.RemoveProperty( "EvoGenFitness");
            }
            if( mol.IsPropertyStored( "EvoGenHistory"))
            {
              mol.RemoveProperty( "EvoGenHistory");
            }
            if( mol.IsPropertyStored( "EvoGenIdentifier"))
            {
              mol.RemoveProperty( "EvoGenIdentifier");
            }
            mol.StoreProperty( "EvoGenFitness", util::Format()( mols[ i].m_Fitness));
            mol.StoreProperty( "EvoGenHistory", util::Format()( mols[ i].m_History));
            mol.StoreProperty( "EvoGenIdentifier", util::Format()( mols[ i].m_Identifier));
            mol.WriteMDL( out);
            out_csv << "0," << mols[ i].m_Fitness << "\n";
          }
          out_csv.flush();
          io::File::CloseClearFStream( out);
        }

        // open the activity log
        std::string activity_log( pop_sdf_prefix + "_activity_log.json");
        optimizer.SetLogFile( activity_log);
        optimizer.StartLogging();

        // execute the optimization until we have iterated enough times
        chemistry::FragmentEnsemble mols_across_generations;
        chemistry::ConstitutionSet unique_mols;
        while( optimizer.GetMolInfos().size() < n_iters)
        {
          // population ID number
          size_t pop_no( optimizer.GetMolInfos().size());

          // check for errors
          if( optimizer.Next() == -1)
          {
            BCL_MessageStd( "Could not optimize population " + util::Format()( optimizer.GetMolInfos().size()) + ", optimizer returned error");
            break;
          }

          // open up an sdf file to write molecule info to
          std::string out_filename( pop_sdf_prefix + "_" + util::Format()( pop_no) + ".sdf.gz");
          const std::vector< MolInfo> &mols( optimizer.GetMolInfos().back());

          // write all molecules, storing necessary properties
          io::File::MustOpenOFStream( out, out_filename);
          for( size_t i( 0), end_i( mols.size()); i < end_i; ++i)
          {
            chemistry::FragmentComplete mol( mols[ i].m_Molecule);
            if( mol.IsPropertyStored( "EvoGenFitness"))
            {
              mol.RemoveProperty( "EvoGenFitness");
            }
            if( mol.IsPropertyStored( "EvoGenHistory"))
            {
              mol.RemoveProperty( "EvoGenHistory");
            }
            if( mol.IsPropertyStored( "EvoGenIdentifier"))
            {
              mol.RemoveProperty( "EvoGenIdentifier");
            }
            mol.StoreProperty( "EvoGenFitness", util::Format()( mols[ i].m_Fitness));
            mol.StoreProperty( "EvoGenHistory", util::Format()( mols[ i].m_History));
            mol.StoreProperty( "EvoGenIdentifier", util::Format()( mols[ i].m_Identifier));
            mol.WriteMDL( out);
            out_csv << pop_no << "," << mols[ i].m_Fitness << "\n";

            // Save current generation to all, but make sure the final selection does not include duplicates carried across generations
            if( unique_mols.Insert( chemistry::FragmentConstitutionShared( mols[ i].m_Molecule)).second)
            {
              mols_across_generations.PushBack( mol);
            }
          }

          // write info to the all_scores csv without closing the stream
          out_csv.flush();

          // close the molecule output file
          io::File::CloseClearFStream( out);

          // display message to user
          BCL_MessageStd
          (
            " Opti has generated " + util::Format()( optimizer.GetMolInfos().size()) + " populations with "
            + util::Format()( optimizer.GetMolInfos().back().size())
          );
        }

        // close down the json log file
        optimizer.StopLogging();
        io::File::CloseClearFStream( out_csv);

        //Sort by fitness and save the top N
        size_t n_best( m_PopulationSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>()), itr_index( 0);
        mols_across_generations.Sort( "EvoGenFitness"); // I cannot figure out how to reverse iterate with a FragmentEnsemble
        storage::List< chemistry::FragmentComplete> sorted_mols( mols_across_generations.GetMolecules());
        io::File::MustOpenOFStream( out, final_output_fn);
        BCL_MessageStd( "Collecting best molecules...");
        for( storage::List< chemistry::FragmentComplete>::reverse_iterator itr( sorted_mols.ReverseBegin()); itr != sorted_mols.ReverseEnd(); ++itr, ++itr_index)
        {
          if( itr_index < n_best)
          {
            itr->WriteMDL( out);
          }
        }

        timer.Stop();
        timer.WriteMessage();

        BCL_MessageStd( "Finished evolving molecules");

        return 0;
      } // Main()

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    private:

      // a static instance of this class
      static const ApplicationType EvoGen_Instance;

    }; // EvoGen

    //! @brief standard constructor
    EvoGen::EvoGen() :
      m_StartingMolsFlag
      (
        new command::FlagStatic
        (
          "starting_mols", "molecules to begin the run",
          command::Parameter
          (
            "filename", "sdf file of starting molecules",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_InsertMolsFlag
      (
        new command::FlagStatic
        (
          "insert_mols", "molecules to insert into the population during the run",
          command::Parameter
          (
            "filename", "sdf file of insertion molecules",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ReactionsFlag
      (
        new command::FlagStatic
        (
          "reactions", "reactions to use for structure modification",
          command::Parameter
          (
            "filename", "rxn file containing desired reactions",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ReagentsFlag
      (
        new command::FlagStatic
        (
          "reagents", "reagents to use in reactions",
          command::Parameter
          (
            "filename", "sdf file of reagnet molecules",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ModelDirectoryFlag
      (
        new command::FlagDynamic
        (
          "model_dir", "QSAR model directory",
          command::Parameter
          (
            "directory", "directory containing models that could be used in a Prediction() descriptor",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ModelPrefixFlag
      (
        new command::FlagDynamic
        (
          "model_prefix", "model prefix",
          command::Parameter
          (
            "prefix", "the prefix for models in the directory specified by -model_dir"
          )
        )
      ),
      m_ExtScoringCommandFlag
      (
        new command::FlagDynamic
        (
          "score_command", "the command to use for external scoring",
          command::Parameter
          (
            "command",
            "the command to use; should return 0 on success, with "
              "usage: program <input sdf> <output sdf>; output sdf should have a field \"EvoGenFitness\" where scores can be read"
          )
        )
      ),
      m_DruglikenessModelFlag
      (
        new command::FlagDynamic
        (
          "druglikeness_descriptor", "a descriptor that provides binary 1/0 output if a molecule is druglike or not",
          command::Parameter
          (
            "descriptor", "A descriptor that returns 1 if a molecule is druglike and 0 otherwise"
          )
        )
      ),
      m_WeightScoreFlag
      (
        new command::FlagDynamic
        (
          "use_weight_term", "whether to use a molecular weight scoring term to restrict molecule size",
          command::Parameter
          (
            "use term", "set this flag to use a weight term"
          )
        )
      ),
      m_OutputDirFlag
      (
        new command::FlagStatic
        (
          "output_dir", "the output directory for the run",
          command::Parameter
          (
            "directory", "the directory for run output",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "output_prefix",
          "what output files should be prefixed with",
          command::Parameter
          (
            "prefix", "identifier for the run",
            "evogen"
          )
        )
      ),
      m_PopulationSizeFlag
      (
        new command::FlagStatic
        (
          "population_size", "size of the population to use",
          command::Parameter
          (
            "size", "number of molecules present in each iteration of the population",
            command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
            "50"
          )
        )
      ),
      m_OversampleFlag
      (
        new command::FlagStatic
        (
          "oversample_factor",
          "the oversampling factor during the run.  the run will generate factor*pop size members during the run before performing replacement",
          command::Parameter
          (
            "factor", "oversampling factor",
            command::ParameterCheckRanged< float>( 1.0, std::numeric_limits< float>::max()),
            "10"
          )
        )
      ),
      m_ShuffleFlag
      (
        new command::FlagDynamic
        (
          "shuffle",
          "whether to randomize initial population",
          command::Parameter
          (
            "shuffle",
            "if set, the initial population will be chosen randomly from the input molecules"
          )
        )
      ),
      m_IterationsFlag
      (
        new command::FlagStatic
        (
          "iterations",
          "the criterion to use to stop the algorithm",
          command::Parameter
          (
            "iterations",
            "the number of iterations to run",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
            "25"
          )
        )
      ),
      m_ActiveCutoffFlag
      (
        new command::FlagDynamic
        (
          "active_cutoff", "the model score value that differentiates actives from inactives",
          command::Parameter
          (
            "cutoff", "the model score cutoff",
            command::ParameterCheckRanged< float>( std::numeric_limits< float>::min(), std::numeric_limits< float>::max())
          )
        )
      ),
      m_SelectionTypeFlag
      (
        new command::FlagDynamic
        (
          "selection_type", "the type of molecule selection to perform",
          command::Parameter
          (
            "selection_type", "selection type",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "tournament", "top")),
            "tournament"
          )
        )
      ),
      m_EvolutionBalanceTypeFlag
      (
        new command::FlagDynamic
        (
          "evolution_balance_type", "the type of balance applied to molecule evolution processes",
          command::Parameter
          (
            "evolution_balance_type", "evolution balance type",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "alchemical_mutate", "reaction_dominant", "reaction_insertion_only", "recombination_dominant", "recombination_insertion_only", "balanced")),
            "reaction_insertion_only"
          )
        )
      ),
      m_ReplaceTournSizeFlag
      (
        new command::FlagDynamic
        (
          "rep_tourn_factor", "the fraction of molecules to consider for replacement tournament selection steps",
          command::Parameter
          (
            "size_factor", "the size factor to use",
            command::ParameterCheckRanged< float>( 0.0, 1.0),
            "0.2"
          )
        )
      ),
      m_ModifyTournSizeFlag
      (
        new command::FlagDynamic
        (
          "mod_tourn_factor", "the fraction of molecules to consider for modification tournament selection steps",
          command::Parameter
          (
            "size_factor", "the size factor to use",
            command::ParameterCheckRanged< float>( 0.0, 1.0),
            "0.2"
          )
        )
      ),
      m_RetirementTypeFlag
      (
        new command::FlagDynamic
        (
          "retirement_type", "the fraction of molecules to consider for modification tournament selection steps",
          command::Parameter
          (
            "type", "the type of retirement to use",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "all", "none", "probabilistic")),
            "none"
          )
        )
      ),
      m_ImplementationFlag
      (
        new command::FlagStatic
        (
          "implementation",
          "method to mutate molecules",
          command::Parameter
          (
            "mutate",
            "",
            command::ParameterCheckSerializable
            (
              util::Implementation< chemistry::FragmentMutateInterface>()
            )
          )
        )
      )
    {
    }

    //! @brief less-than operator for ShPtrs to MolInfos
    //! @return true if left operand is less than right operand
    bool operator <( const util::ShPtr< MolInfo> &FIRST, const util::ShPtr< MolInfo> &SECOND)
    {
      BCL_Assert( FIRST.IsDefined() && SECOND.IsDefined(), "comparison would have dereferenced null pointer");
      return *FIRST < *SECOND;
    }

    //! @brief greater-than operator for ShPtrs to MolInfos
    //! @return true if left operand is greater than right operand
    bool operator >( const util::ShPtr< MolInfo> &FIRST, const util::ShPtr< MolInfo> &SECOND)
    {
      BCL_Assert( FIRST.IsDefined() && SECOND.IsDefined(), "comparison would have dereferenced null pointer");
      return *FIRST > *SECOND;
    }

    //! @brief less-than operator for pair of MolInfo* and size_t
    //! @return true if left operand is less than right operand
    bool operator <( const std::pair< const MolInfo *, size_t> &FIRST, const std::pair< const MolInfo *, size_t> &SECOND)
    {
      return FIRST.first && SECOND.first ? *FIRST.first < *SECOND.first : 0;
    }

    //! @brief greater-than operator for pair of MolInfo* and size_t
    //! @return true if left operand is greater than right operand
    bool operator >( const std::pair< const MolInfo *, size_t> &FIRST, const std::pair< const MolInfo *, size_t> &SECOND)
    {
      return FIRST.first && SECOND.first ? *FIRST.first > *SECOND.first : 0;
    }

    // static instance of this class used for adding to the command line
    const ApplicationType EvoGen::EvoGen_Instance
    (
      GetAppGroups().AddAppToGroup( new EvoGen(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app
} // namespace bcl
