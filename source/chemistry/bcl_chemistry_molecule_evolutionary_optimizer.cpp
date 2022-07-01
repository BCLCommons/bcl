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
#include <chemistry/bcl_chemistry_fragment_split_gadd_fragments.h>
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// unit header
#include "chemistry/bcl_chemistry_molecule_evolutionary_optimizer.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_configuration_set.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_evolve_implementations.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_graph_marker.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_subgraph.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "math/bcl_math_running_average.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_stopwatch.h"

// external includes

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeEvolutionaryOptimizer::MoleculeEvolutionaryOptimizer() :
      m_GenerateMax( 10),
      m_FinalPopSize( 10),
      m_MaxFailedAttempts( 10),
      m_Operations( 3, size_t( 0)),
      m_SelectionType( e_Top),
      m_ModelType( e_Internal),
      m_ModelCmd( std::string()),
      m_RetirementType( e_DoNotRetire),
      m_EvolutionBalanceType( e_ReactionDominant),
      m_RecombineOp( FragmentEvolveImplementations::EvolveType::e_Combine, false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoleculeEvolutionaryOptimizer
    MoleculeEvolutionaryOptimizer *MoleculeEvolutionaryOptimizer::Clone() const
    {
      return new MoleculeEvolutionaryOptimizer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoleculeEvolutionaryOptimizer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of this class
    //! @return the name of this class
    const std::string &MoleculeEvolutionaryOptimizer::GetAlias() const
    {
      static const std::string s_name( "MoleculeEvolutionaryOptimizer");
      return s_name;
    }

    //! @brief get the internal molecule data
    //! @return a vector containing population data for each iteration so far
    const std::vector< std::vector< MoleculeEvolutionInfo> > &MoleculeEvolutionaryOptimizer::GetMoleculeEvolutionInfos() const
    {
      return m_Populations;
    }

    //! @brief Get filename for EvoGen log file
    const std::string &MoleculeEvolutionaryOptimizer::GetLogFile() const
    {
      return m_LogFile;
    }

    //! @brief Get the final size of each population.
    const size_t MoleculeEvolutionaryOptimizer::GetFinalPopSize() const
    {
      return m_FinalPopSize;
    }

    //! @brief Get the maximum number of molecules to generate during each iteration.
    const size_t MoleculeEvolutionaryOptimizer::GetMaxToGenerate() const
    {
      return m_GenerateMax;
    }

    //! @brief Get molecule selection type to keep highest-scoring molecules
    const MoleculeEvolutionaryOptimizer::SelectionType &MoleculeEvolutionaryOptimizer::GetSelectionType() const
    {
      return m_SelectionType;
    }

    //! @brief Get molecule replacement method to use tournament selection
    const float MoleculeEvolutionaryOptimizer::GetReplacementTypeTournament() const
    {
      return m_ReplacementTournSizeFactor;
    }

    //! @brief Get molecule modification method to use tournament selection
    const float MoleculeEvolutionaryOptimizer::GetModifyTypeTournament() const
    {
      return m_ModifyTournSizeFactor;
    }

    //! @brief Get retirement policy so that all parents are discarded
    const MoleculeEvolutionaryOptimizer::ParentRetirementType &MoleculeEvolutionaryOptimizer::GetRetirementType() const
    {
      return m_RetirementType;
    }

    //! @brief Get molecule evolution type to primarily perform balanced processes
    const MoleculeEvolutionaryOptimizer::EvolutionBalance &MoleculeEvolutionaryOptimizer::GetEvolutionBalanceType() const
    {
      return m_EvolutionBalanceType;
    }

    //! @brief Use BCL-internal models for molecule scoring
    const MoleculeEvolutionaryOptimizer::ModelType &MoleculeEvolutionaryOptimizer::GetModelType() const
    {
      return m_ModelType;
    }

    //! @brief Get the reaction operation structure
    const std::string &MoleculeEvolutionaryOptimizer::GetReactionOperationLabel() const
    {
      return m_ReactOp.GetAlias();
    }

    //! @brief Get the one-shot reaction operation structure
    const std::string &MoleculeEvolutionaryOptimizer::GetAlchemicalOperationLabel() const
    {
      return m_Mutate.GetAlias();
    }

    //! @brief Get the descriptor to use for scoring, if internal models are used
    const std::string &MoleculeEvolutionaryOptimizer::GetModelDescriptorLabel() const
    {
      return m_Scorer.GetAlias();
    }

    //! @brief Get the descriptor to use for scoring, if internal models are used
    const descriptor::CheminfoProperty &MoleculeEvolutionaryOptimizer::GetModelDescriptor() const
    {
      return m_Scorer;
    }

    //! @brief Get the external script path if external scoring is used
    const std::string &MoleculeEvolutionaryOptimizer::GetModelCmd()
    {
      return m_ModelCmd;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set filename for EvoGen log file
    //! @details this file is a json-formatted log file containing information about the generated molecules
    void MoleculeEvolutionaryOptimizer::SetLogFile( const std::string &FILENAME)
    {
      m_LogFile = FILENAME;
    }

    //! @brief set up the initial population from an ensemble of molecules
    //! @details molecules must have >= 1 atoms or they will be discarded.
    //! @param MOLS the structures to use, will be copied
    void MoleculeEvolutionaryOptimizer::SetInitialPopulation( const FragmentEnsemble &MOLS)
    {
      // check necessary components are available
      BCL_Assert( m_Scorer.IsDefined() || !m_ModelCmd.empty(), "Must set score descriptor before calling SetInitialPopulation");

      // clear all data and start with a fresh, empty, population
      m_Populations.clear();
      m_Populations.push_back( std::vector< MoleculeEvolutionInfo>());
      std::vector< MoleculeEvolutionInfo> &new_pop( m_Populations[ 0]);

      // perform checks on each molecule and add to new population with proper information
      size_t mol_no( 0);
      for
      (
        FragmentEnsemble::const_iterator itr_mol( MOLS.Begin()), itr_mol_end( MOLS.End());
        itr_mol != itr_mol_end;
        ++itr_mol
      )
      {
        if( itr_mol->GetNumberAtoms())
        {
          new_pop.push_back
          (
            MoleculeEvolutionInfo( "P" + util::Format()( 0) + "M" + util::Format()( mol_no++), *itr_mol) //, 0, "Begin")
          );
        }
      }

      ScoreMolecules( new_pop);

      // sort molecules by fitness low->high
      std::sort( new_pop.begin(), new_pop.end());
    }

    //! @brief set the final size of each population.
    //! @param POP_SIZE the final size of populations
    void MoleculeEvolutionaryOptimizer::SetFinalPopSize( const size_t &POP_SIZE)
    {
      m_FinalPopSize = POP_SIZE;
    }

    //! @brief set the maximum number of molecules to generate during each iteration.  This will be pruned to
    //!  m_FinalPopSize afterwards
    //! @param MAX maximum number of molecules to generate
    void MoleculeEvolutionaryOptimizer::SetMaxToGenerate( const size_t &MAX)
    {
      m_GenerateMax = MAX;
    }

    //! @brief set molecule selection type to keep highest-scoring molecules
    void MoleculeEvolutionaryOptimizer::SetSelectionTop()
    {
      m_SelectionType = e_Top;
      BCL_MessageStd( "Selection type is: " + util::Format()( m_SelectionType));
    }

    //! @brief Set molecule replacement method to use tournament selection
    //! @param FACTOR the percentage of the available data that will be used in a single tournament round
    //! @details will warn if FACTOR < 0 or FACTOR > 1 and set FACTOR to 0.0 or 1.0, respectively.  The number of
    //!  molecules selected will be 1 <= NUMBER <= data size regardless of FACTOR
    void MoleculeEvolutionaryOptimizer::SetReplacementTypeTournament( const float &FACTOR)
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
    void MoleculeEvolutionaryOptimizer::SetModifyTypeTournament( const float &FACTOR)
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
    void MoleculeEvolutionaryOptimizer::SetEvolutionBalanceAlchemicalMutate()
    {
      m_EvolutionBalanceType = e_AlchemicalMutate;
      BCL_MessageStd( "Evolution balance type is: AlchemicalMutate");
    }

    void MoleculeEvolutionaryOptimizer::SetEvolutionBalanceReactionDominant()
    {
      m_EvolutionBalanceType = e_ReactionDominant;
      BCL_MessageStd( "Evolution balance type is: ReactionDominant");
    }

    //! @brief set molecule evolution type to primarily perform reactions and insertions
    void MoleculeEvolutionaryOptimizer::SetEvolutionBalanceReactionInsertionOnly()
    {
      m_EvolutionBalanceType = e_ReactionInsertionOnly;
      BCL_MessageStd( "Evolution balance type is: ReactionInsertionOnly");
    }

    //! @brief set molecule evolution type to primarily perform recombinations
    void MoleculeEvolutionaryOptimizer::SetEvolutionBalanceRecombinationDominant()
    {
      m_EvolutionBalanceType = e_RecombinationDominant;
      BCL_MessageStd( "Evolution balance type is: RecombinationDominant");
    }

    //! @brief set molecule evolution type to primarily perform recombinations and insertions
    void MoleculeEvolutionaryOptimizer::SetEvolutionBalanceRecombinationInsertionOnly()
    {
      m_EvolutionBalanceType = e_RecombinationInsertionOnly;
      BCL_MessageStd( "Evolution balance type is: RecombinationInsertionOnly");
    }

    //! @brief set molecule evolution type to primarily perform balanced processes
    void MoleculeEvolutionaryOptimizer::SetEvolutionBalancedBalanced()
    {
      m_EvolutionBalanceType = e_ReactionRecombinationBalanced;
      BCL_MessageStd( "Evolution balance type is: ReactionRecombinationBalanced");
    }

    //! @brief use BCL-internal models for molecule scoring
    void MoleculeEvolutionaryOptimizer::SetModelTypeInternal()
    {
      m_ModelType = e_Internal;
    }

    //! @brief use external script for molecule scoring
    void MoleculeEvolutionaryOptimizer::SetModelTypeExternal()
    {
      m_ModelType = e_External;
    }

    //! @brief initialize the reaction operation structure
    //! @brief REACTANT_FILENAME filename from which to read reactant molecules
    //! @brief REACTION_DIRNAME directory in which RXN files should be found
    void MoleculeEvolutionaryOptimizer::SetupReactOperation( const std::string &REACTANT_FILENAME, const std::string &REACTION_DIRNAME)
    {
      m_ReactOp = FragmentReact( ReactionSearch( REACTANT_FILENAME, REACTION_DIRNAME));
    }

    //! @brief initialize the one-shot reaction operation structure
    //! @brief IMPLEMENTATION alchemical mutate implementation
    void MoleculeEvolutionaryOptimizer::SetupAlchemicalMutate( const std::string &IMPLEMENTATION)
    {
      // load fragment pool
      m_Mutate = IMPLEMENTATION;
    }

    //! @brief set the descriptor to use for scoring, if internal models are used
    //! @param DESCRIPTOR a string which should property encode a descriptor::CheminfoProperty
    //! @details this asserts that internal models are used to prevent programming errors
    void MoleculeEvolutionaryOptimizer::SetModelDescriptor( const std::string &DESCRIPTOR)
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
    void MoleculeEvolutionaryOptimizer::SetModelCmd( const std::string &CMD)
    {
      BCL_Assert
      (
        m_ModelType == e_External,
        "Programming error, tried setting a model script when using BCL-internal models!"
      );
      m_ModelCmd = CMD;
    }

    //! @brief set retirement policy so that no parents are discarded
    void MoleculeEvolutionaryOptimizer::SetRetirementTypeNone()
    {
      m_RetirementType = e_DoNotRetire;
    }

    //! @brief set retirement policy so that parents are discarded probabilistically based on age
    void MoleculeEvolutionaryOptimizer::SetRetirementTypeProbabilistic()
    {
      m_RetirementType = e_RetireProbabilisticByAge;
    }

    //! @brief set retirement policy so that all parents are discarded
    void MoleculeEvolutionaryOptimizer::SetRetirementTypeAll()
    {
      m_RetirementType = e_RetireAll;
    }

    //! @brief set up the molecule insertion object
    //! @param INSERT_MOL_FILENAME sdf files which should be added to populations during a run
    void MoleculeEvolutionaryOptimizer::SetupInsertOperation( const std::string &INSERT_MOL_FILENAME)
    {
      io::IFStream in;
      io::File::MustOpenIFStream( in, INSERT_MOL_FILENAME);
      FragmentFeed ff( in, sdf::e_Saturate);
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
    size_t MoleculeEvolutionaryOptimizer::SelectMoleculeTournament
    (
      const std::vector< MoleculeEvolutionInfo> &MOLS,
      size_t TOURN_SIZE,
      const std::set< int> &IGNORE_INDICES
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

      // make a vector of MoleculeEvolutionInfo* and save the associated index in MOLS so that we can sort them by fitness
      // track fitness sums and minimum to properly calculate probabilities
      std::vector< std::pair< const MoleculeEvolutionInfo *, size_t> > m_ptrs;
      double total_fitness( 0);
      for( size_t i( 0); i < inds.size(); ++i)
      {
        m_ptrs.push_back( std::pair< const MoleculeEvolutionInfo *, size_t>( &MOLS[ inds[ i]], inds[ i]));
        total_fitness += MOLS[ inds[ i]].GetMoleculeFitness();
      }

      // sort the pointers from high to low by fitness
      std::sort( m_ptrs.begin(), m_ptrs.end(), std::greater< std::pair< const MoleculeEvolutionInfo *, size_t> >());

      // make sure total fitness is not too small for division purposes
      total_fitness = std::max< double>( 0.001, total_fitness);

      // probability of selecting the first member is proportional to its fitness relative to the whole data set,
      // but should not be lower than 10%
      double p( std::max< double>( 0.1, m_ptrs[ 0].first->GetMoleculeFitness() / total_fitness));
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

    //! @brief helper funciton to score a vector of MoleculeEvolutionInfos
    //! @param MOL_INFOS the MoleculeEvolutionInfos containing molecules to score
    //! @details this updates the MoleculeEvolutionInfo.m_Fitness function of all molecules
    void MoleculeEvolutionaryOptimizer::ScoreMolecules( std::vector< MoleculeEvolutionInfo> &MOL_INFOS)
    {
      switch( m_ModelType)
      {
        case e_Internal:
          for
          (
            std::vector< MoleculeEvolutionInfo>::iterator itr_mi( MOL_INFOS.begin()), itr_mi_end( MOL_INFOS.end());
            itr_mi != itr_mi_end;
            ++itr_mi
          )
          {
            itr_mi->SetMoleculeFitness( ScoreMoleculeInternal( itr_mi->GetMolecule()));
          }
          break;
        case e_External:
          ScoreMoleculesExternal( MOL_INFOS);
          break;
        default:
          BCL_Exit( "Neither internal nor external model was specified", -1);
      }
    }

    //! @brief helper function for scoring a single molecule using internal BCL models
    //! @param MOL the molecule to score
    //! @return the score of the molecule
    //! @details does not do any preprocessing of the molecule, this is up to the caller
    float MoleculeEvolutionaryOptimizer::ScoreMoleculeInternal( const FragmentComplete &MOL) const
    {
      BCL_Assert( m_Scorer.IsDefined(), "scorer was not defined");
      return m_Scorer->SumOverObject( MOL)( 0);
    }

    //! @brief helper function for scoring a single molecule using an external script
    //! @param MOL_INFOS the molecules to score
    //! @details this will update the MOL_INFOS structure directly with the new fitness scores as read from the
    //!  output of the called script
    void MoleculeEvolutionaryOptimizer::ScoreMoleculesExternal( std::vector< MoleculeEvolutionInfo> &MOL_INFOS) const
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
      std::map< std::string, MoleculeEvolutionInfo *> mol_names;

      // write out molecule data to sdf
      io::OFStream out;
      io::File::MustOpenOFStream( out, filename_out);
      size_t mol_no( 0);
      for
      (
        std::vector< MoleculeEvolutionInfo>::iterator itr_mi( MOL_INFOS.begin()), itr_mi_end( MOL_INFOS.end());
        itr_mi != itr_mi_end;
        ++itr_mi, ++mol_no
      )
      {
        BCL_Assert
        (
          !itr_mi->GetMoleculeIdentifier().empty(),
          "Programming error: molecule " + util::Format()( mol_no) + " has an empty identifier!"
        );

        // remove stale properties
        if( !itr_mi->GetMolecule().IsPropertyStored( "TmpEvoGenIdentifier"))
        {
          itr_mi->GetMoleculeNonConst().RemoveProperty( "TmpEvoGenIdentifier");
        }

        // add an identifier to the molecule
        itr_mi->GetMoleculeNonConst().StoreProperty( "TmpEvoGenIdentifier", itr_mi->GetMoleculeIdentifier());

        // All identifiers should be unique so this should never actually happen, but double check
        if( mol_names.find( itr_mi->GetMoleculeIdentifier()) != mol_names.end())
        {
          BCL_MessageStd( "Multiple molecules were given with ID \"" + itr_mi->GetMoleculeIdentifier() + "\"!");
        }
        mol_names[ itr_mi->GetMoleculeIdentifier()] = &( *itr_mi);

        itr_mi->GetMolecule().WriteMDL( out);
      }

      io::File::CloseClearFStream( out);

      // execute the external command: cmd <infile to prog> <outfile from prog>
      std::string cmd( m_ModelCmd + " " + filename_out + " " + filename_in);
      int err = system( cmd.c_str());
      BCL_Assert( err == 0, "Could not execute command \"" + cmd + "\", exitted with status " + util::Format()( err));

      // read in what cmd wrote out
      io::IFStream in;
      io::File::MustOpenIFStream( in, filename_in);
      FragmentEnsemble in_ens( in, sdf::HydrogenHandlingPref::e_Saturate);

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
      std::map< std::string, MoleculeEvolutionInfo *>::iterator map_end( mol_names.end());
      for
      (
        FragmentEnsemble::const_iterator itr_mol( in_ens.Begin()), itr_mol_end( in_ens.End());
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
        std::map< std::string, MoleculeEvolutionInfo *>::iterator itr_map( mol_names.find( m_name));
        BCL_Assert
        (
          itr_map != map_end,
          "Could not find a molecule with identifier \"" + m_name + "\" in the MoleculeEvolutionInfo vector (out of " + util::Format()( mol_names.size()) + ")"
        );

        // update the molecule with whatever was provided by the scoring script
        itr_map->second->SetMolecule( *itr_mol);
        itr_map->second->SetMoleculeFitness( score);
        itr_map->second->GetMoleculeNonConst().RemoveProperty( "TmpEvoGenIdentifier");
      }
      remove( filename_out.c_str());
      remove( filename_in.c_str());
    }

    //! @brief copy the highest scoring molecules from one population to a number
    //! @param FROM the population to copy from
    //! @param TO the population to copy to
    //! @param NUM the number to copy
    void MoleculeEvolutionaryOptimizer::CopyHighestScoring( std::vector< MoleculeEvolutionInfo> &FROM, std::vector< MoleculeEvolutionInfo> &TO, size_t NUM)
    {
      std::sort( FROM.begin(), FROM.end());
      size_t n_added( 0);
      for
      (
        std::vector< MoleculeEvolutionInfo>::const_reverse_iterator itr( FROM.rbegin()), itr_end( FROM.rend());
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
    void MoleculeEvolutionaryOptimizer::CopyByTournament
    (
      const std::vector< MoleculeEvolutionInfo> &FROM,
      std::vector< MoleculeEvolutionInfo> &TO,
      size_t NUM,
      size_t TOURN_SIZE,
      const bool &REPLACEMENT
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

    //! @brief evaluates a descriptor to determine whether a molecule
    //! passes or fails the druglikeness filter
    bool MoleculeEvolutionaryOptimizer::EvaluateDruglikeness( const MoleculeEvolutionInfo &MOL) const
    {
      // evaluate the properties across the molecule
      const linal::Vector< float> lhs_property( m_DruglikenessFilter.First()->SumOverObject( MOL.GetMolecule()));
      const linal::Vector< float> rhs_property( m_DruglikenessFilter.Third()->SumOverObject( MOL.GetMolecule()));

      // if no value is returned then assume
      if( !lhs_property.GetSize() || !rhs_property.GetSize())
      {
        return false;
      }

      // perform comparison
      return ( **m_DruglikenessFilter.Second())( lhs_property( 0), rhs_property( 0));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief execute a reaction/addition operation to generate new molecules
    //! @brief MOLS the vector to add the new molecules to
    //! @details this only generates structures and sets history values appropriately, but does nothing for ensuring
    //!  molecules are suitable for scoring or that their
    void MoleculeEvolutionaryOptimizer::GenerateMolecules( std::vector< MoleculeEvolutionInfo> &MOLS) const
    {
      if( m_Populations.empty())
      {
        return;
      }

      // reference to most recent population, should be sorted by score l->h
      const std::vector< MoleculeEvolutionInfo> &last_pop( m_Populations.back());
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
//        const FragmentComplete &picked_mol( last_pop[ mol_no].GetMolecule());
        FragmentComplete picked_mol( last_pop[ mol_no].GetMolecule());

        // mutate
        auto product( ( *m_Mutate)( picked_mol));

        // update info
        if( product.GetArgument().IsDefined())
        {
          auto hist( last_pop[ mol_no].GetMoleculeHistory());
          hist.Append( "AlchemicalMutate," + m_Mutate->GetAlias() + last_pop[ mol_no].GetMoleculeIdentifier());
          MOLS.push_back( MoleculeEvolutionInfo( "", *( product.GetArgument()), util::GetUndefined< float>(), hist));
        }
      }

      // rxn-type reaction
      else if( op_test < react)
      {
        ++m_Operations[ 0];

        // select molecule via tournament selection
        size_t mol_no( SelectMoleculeTournament( last_pop, tourn_size));

        // react the molecule with something random
        const FragmentComplete &picked_mol( last_pop[ mol_no].GetMolecule());
        storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> res
        (
          m_ReactOp.ReactRandom( picked_mol)
        );

        // a reaction may have generated more than one molecule, add each of them to the new population
        for
        (
          FragmentEnsemble::const_iterator itr_mol( res.Second().Begin()), itr_mol_end( res.Second().End());
            itr_mol != itr_mol_end;
            ++itr_mol
        )
        {
//          std::stringstream hist;
          auto hist( last_pop[ mol_no].GetMoleculeHistory());
          hist.Append( "React," + RemoveWhitespace( res.First()->GetDescription()) + "," + last_pop[ mol_no].GetMoleculeIdentifier());
//          hist << "React,";
//          hist << RemoveWhitespace( res.First()->GetDescription());
//          hist << "," << last_pop[ mol_no].GetMoleculeIdentifier();
          MOLS.push_back( MoleculeEvolutionInfo( "", *itr_mol, util::GetUndefined< float>(), hist));
        }
      }
      // recombination
      else if( op_test < recombine)
      {
        ++m_Operations[ 1];

        // select molecule from tournament selection
        size_t mol_one( SelectMoleculeTournament( last_pop, tourn_size));
        const FragmentComplete &picked_mol_one( last_pop[ mol_one].GetMolecule());
        util::ShPtrVector< FragmentComplete> new_mols;

        // select recombination type randomly
        size_t recomb_type( random::GetGlobalRandom().Random< size_t>( 0, 2));
        if( recomb_type == 0 || recomb_type == 1)
        {
          // combine selected molecule with random molecule from the inserted (or "migrant") mols population
          m_RecombineOp = FragmentEvolveImplementations( FragmentEvolveImplementations::e_Combine, false);
          size_t mol_two( random::GetGlobalRandom().Random< size_t>( 0, m_InsertMols.GetSize() - 1));
          const FragmentComplete &picked_mol_two( *m_InsertMols( mol_two));
          new_mols = m_RecombineOp.Combine( picked_mol_one, picked_mol_two, 50);
        }
        else if( recomb_type == 3)
        {
          // split GADD fragments from a random inserted (or "migrant) molecule and
          m_RecombineOp = FragmentEvolveImplementations( FragmentEvolveImplementations::e_FragAdd, false);
//              size_t mol_two( random::GetGlobalRandom().Random< size_t>( 0, m_InsertMols.GetSize() - 1));
//              const FragmentComplete &picked_mol_two( *m_InsertMols( mol_two));

          //List of component fragment indices
          storage::List< storage::Vector< size_t> > gadd_fragment_indices;

          // GADD fragment rings
          if( random::GetGlobalRandom().Random< size_t>( 0, 1) == 0)
          {
            FragmentSplitGADDFragments splitter;
//                splitter.GetComponentVertices( picked_mol_two, picked_mol_two.)
          }

          // GADD fragment chains
          else
          {
            FragmentSplitGADDFragments splitter;
          }

//              new_mols = m_RecombineOp.MutateAdd()

        }
        else if( recomb_type == 2)
        {
          m_RecombineOp = FragmentEvolveImplementations( FragmentEvolveImplementations::e_FragDel, false);
          new_mols = m_RecombineOp.MutateDel( picked_mol_one);
        }

        // now add the new compounds to the mix
        for
        (
          util::ShPtrVector< FragmentComplete>::iterator new_mols_itr( new_mols.Begin()), new_mols_itr_end( new_mols.End());
            new_mols_itr != new_mols_itr_end;
            ++new_mols_itr
        )
        {
//          std::stringstream hist;
          auto hist( last_pop[mol_one].GetMoleculeHistory());
          hist.Append( "Recombine,");
//          hist << "Recombine,";
//          hist << last_pop[ mol_one].GetMoleculeHistory();
          MOLS.push_back( MoleculeEvolutionInfo( "", **new_mols_itr, util::GetUndefined< float>(), hist));
        }
      }

      else
      {
        ++m_Operations[ 2];
        if( m_InsertMols.GetSize() > 1)
        {
          size_t mol_no( random::GetGlobalRandom().Random< size_t>( 0, m_InsertMols.GetSize() - 1));
          MOLS.push_back( MoleculeEvolutionInfo());
          MoleculeEvolutionInfo &mol = MOLS.back();
          mol.SetMolecule( *m_InsertMols( mol_no));
          mol.AppendToMoleculeHistory( "Insert,M" + util::Format()( mol_no));
        }
        else if( m_InsertMols.GetSize() == 1)
        {
          MOLS.push_back( MoleculeEvolutionInfo());
          MoleculeEvolutionInfo &mol = MOLS.back();
          mol.SetMolecule( *m_InsertMols( 0));
          mol.AppendToMoleculeHistory( "Insert,M" + util::Format()( 0));
        } // else noop
      }
    }

    //! @brief execute an iteration of molecule generation/scoring/replacement
    //! @return 0 on success, negative value on error
    int MoleculeEvolutionaryOptimizer::Next()
    {
      std::vector< size_t> prev_ops( m_Operations);
      std::vector< MoleculeEvolutionInfo> next_pop;
      next_pop.reserve( m_GenerateMax + ( m_Populations.empty() ? 0 : m_Populations.back().size()));
      ConstitutionSet unique_mols;

      size_t tries = 0;
      size_t mol_no = 0;

      // do this until we reach the specified number of molecules
      while( next_pop.size() < m_GenerateMax)
      {

        std::vector< MoleculeEvolutionInfo> mols;
        GenerateMolecules( mols);

        // add the newly generated molecules to the population one by one...
        for( size_t i( 0); i < mols.size(); ++i)
        {

          util::ShPtr< FragmentComplete> final_mol;
          m_Mutate.IsDefined() ?
          final_mol = util::ShPtr< FragmentComplete> ( util::ShPtr< FragmentComplete>( new FragmentComplete( mols[i].GetMolecule()))) :
          final_mol = util::ShPtr< FragmentComplete> ( FragmentEvolveBase::FinalizeMolecule( mols[ i].GetMolecule()));
          if( final_mol.IsDefined())
          {
            mols[ i].SetMolecule( *final_mol);

            // make sure we don't generate duplicates in a single population
            if
            (
                mols[ i].GetMolecule().GetNumberAtoms() &&
                unique_mols.Insert( FragmentConstitutionShared( mols[ i].GetMolecule())).second &&
                EvaluateDruglikeness( mols[ i])
            )
            {
              mols[ i].SetMoleculeIdentifier( "P" + util::Format()( m_Populations.size()) + "M" + util::Format()( mol_no++));
              next_pop.push_back( mols[ i]);
              tries = 0;
            }
          }
        }

        // if we are failing really hard, then we should probably just abort...
        BCL_Assert( ++tries < m_MaxFailedAttempts, "Tried to generate a unique molecule, took more than 100 tries, bailing out");

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

      storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double>> refined_mols( next_pop.size());

      // score all molecules as a batch
      BCL_MessageStd( "Scoring molecules");
      ScoreMolecules( next_pop);

//      // Write information to the log file in json format.  This is ugly as sin and should be re-done in the future
//      WriteLog( "{ \"pop_number\": " + util::Format()( m_Populations.size()) + ",\n");
//      WriteLog( "  \"pre_replacement\": [\n");
//      for( size_t i( 0); i < next_pop.size(); ++i)
//      {
//        WriteLog( "      { \"id\": \"" + next_pop[ i].GetMoleculeIdentifier() + "\",\n");
//        WriteLog( "        \"history\": \"" + EscapeQuotes( next_pop[ i].GetMoleculeHistory()) + "\",\n");
//        WriteLog( "        \"age\": " + util::Format()( next_pop[ i].GetMoleculeAge()) + ",\n");
//        WriteLog( "        \"score\": " + util::Format()( next_pop[ i].GetMoleculeFitness()) + "}");
//        if( i < ( next_pop.size() - 1))
//        {
//          WriteLog( ",\n");
//        }
//      }
//      WriteLog( "  ],\n");

      // parent replacement: add members from the old population to the new population
      if( m_RetirementType != e_RetireAll)
      {
        BCL_MessageStd( "Copying parents to the new population");
//        WriteLog( "  \"copied_parents\": [\n");

        const std::vector< MoleculeEvolutionInfo> &parent_pop( m_Populations.back());
        for( size_t i( 0); i < parent_pop.size(); ++i)
        {
          double r( 0.5); // test value for keeping a molecule; initial value is a dummy value for Retire None policy
          double cutoff( 1.0); // test threshold; r must be below this value to keep the molecule

          // if retirement policy is probabilistic then "roll the dice" to see if this molecule survives
          if( m_RetirementType == e_RetireProbabilisticByAge)
          {
            r = random::GetGlobalRandom().Random< double>( 0, 1.0); // assign a real test value
            cutoff = exp( double( -0.5) * parent_pop[ i].GetMoleculeAge()); // assign a real test threshold value based on age
          }

          // make sure the molecule passed the test, and it's a unique structure
          if( r < cutoff && unique_mols.Insert( FragmentConstitutionShared( parent_pop[ i].GetMolecule())).second)
          {
            //BCL_MessageStd( "  Keeping parent " + parent_pop[ i].GetMoleculeIdentifier());
            next_pop.push_back( parent_pop[ i]);
            next_pop.back().IncrementMoleculeAge();

//            // write to the log file
//            WriteLog( "      { \"id\": \"" + parent_pop[ i].GetMoleculeIdentifier() + "\",\n");
//            WriteLog( "        \"history\": \"" + EscapeQuotes( parent_pop[ i].GetMoleculeHistory()) + "\",\n");
//            WriteLog( "        \"age\": " + util::Format()( parent_pop[ i].GetMoleculeAge()) + ",\n");
//            WriteLog( "        \"score\": " + util::Format()( parent_pop[ i].GetMoleculeFitness()) + "}");
//            if( i < ( parent_pop.size() - 1))
//            {
//              WriteLog( ",\n");
//            }

          }
          else
          {
            // status message to make users happy
            BCL_MessageStd( "  Retiring parent " + parent_pop[ i].GetMoleculeIdentifier());
          }
        }
        // additional closing braces to make sure json format is ok
//        WriteLog( "\n  ],\n");
      }

      size_t new_pop_no( m_Populations.size()); // the new population number (same as current size since we add one)
      m_Populations.push_back( std::vector< MoleculeEvolutionInfo>()); // allocate space for new population
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

//      // write more log information...
//      WriteLog( "  \"post_replacement\": [\n");
//      for( size_t i( 0); i < m_Populations[ new_pop_no].size(); ++i)
//      {
//        WriteLog( "      { \"id\": \"" + m_Populations[ new_pop_no][ i].GetMoleculeIdentifier() + "\",\n");
//        WriteLog( "        \"history\": \"" + EscapeQuotes( m_Populations[ new_pop_no][ i].GetMoleculeHistory()) + "\",\n");
//        WriteLog( "        \"age\": " + util::Format()( m_Populations[ new_pop_no][ i].GetMoleculeAge()) + ",\n");
//        WriteLog( "        \"score\": " + util::Format()( m_Populations[ new_pop_no][ i].GetMoleculeFitness()) + "}");
//        if( i < ( m_Populations[ new_pop_no].size() - 1))
//        {
//          WriteLog( ",\n");
//        }
//      }
//      WriteLog( "\n  ]\n}\n");

      // tell the user how many reactions/insertions were done
      BCL_MessageStd( "Ran " + util::Format()( m_Operations[ 0] - prev_ops[ 0]) + " reactions");
      BCL_MessageStd( "Ran " + util::Format()( m_Operations[ 1] - prev_ops[ 1]) + " recombinations");
      BCL_MessageStd( "Ran " + util::Format()( m_Operations[ 2] - prev_ops[ 2]) + " inserts");
      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

//    //! @brief open log file for writing; continues if file cannot be opened
//    void MoleculeEvolutionaryOptimizer::StartLogging()
//    {
//      if( !io::File::TryOpenOFStream( m_LogStream, m_LogFile))
//      {
//        BCL_MessageStd( "Could not open " + m_LogFile + " for writing");
//      }
//    }
//
//    //! @brief close/flush logging file stream
//    void MoleculeEvolutionaryOptimizer::StopLogging()
//    {
//      if( m_LogStream.is_open())
//      {
//        io::File::CloseClearFStream( m_LogStream);
//      }
//    }

    //! @brief remove whitespace (via isspace) from a string
    //! @param STR the string to remove whitespace from
    //! @return STR without any whitespace
    std::string MoleculeEvolutionaryOptimizer::RemoveWhitespace( const std::string &STR) const
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
    std::string MoleculeEvolutionaryOptimizer::PrepareForCSV( const std::string &STR) const
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
    std::string MoleculeEvolutionaryOptimizer::EscapeQuotes( const std::string &STR) const
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

//    //! @brief write data to the json log file
//    //! @param STR the string to write
//    void MoleculeEvolutionaryOptimizer::WriteLog( const std::string &STR)
//    {
//      if( m_LogStream.is_open())
//      {
//        m_LogStream << STR;
//      }
//    }

    //! @brief Set the members with LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeEvolutionaryOptimizer::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return true;
    }

    io::Serializer MoleculeEvolutionaryOptimizer::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "Evolves a molecule.");

      // TODO: is it easier for people if this is one string or if i split it into 3?
      member_data.AddOptionalInitializer
      (
        "druglikeness_filter",
        "comparison to determine whether a molecule is druglike; "
        "formatted as a comparison: <lhs_property> <comparison> <rhs_property>",
        io::Serialization::GetAgent( &m_DruglikenessFilterStr)
      );

      return member_data;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoleculeEvolutionaryOptimizer::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &MoleculeEvolutionaryOptimizer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
