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
#include "bcl_app_build_scaffold_library.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_csi_substructure.h"
#include "io/bcl_io_file.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {

    const ApplicationType BuildScaffoldLibrary::BuildScaffoldLibrary_Instance
    (
      GetAppGroups().AddAppToGroup( new BuildScaffoldLibrary(), GetAppGroups().e_Molecule)
    );

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! default constructor
    //! @brief standard constructor
    BuildScaffoldLibrary::BuildScaffoldLibrary() :
        m_EnsembleFileFlag
        (
          new command::FlagStatic
          (
            "ensemble",
            "sdf file containing the ensembles from which to find scaffolds",
            command::Parameter
            (
              "input",
              "an sdf input file",
              command::ParameterCheckFileExistence()
            )
          )
        ),
        m_SamplingFractionFlag
        (
          new command::FlagStatic
          (
            "sampling_fraction",
            "sampling_fraction: fraction of all pairs of molecules that are tested for a common scaffold (max 1)",
            command::Parameter
            (
              "sampling_fraction",
              "fraction of all pairs of molecules that are tested for a common scaffold (max 1)",
              command::ParameterCheckRanged< double>( 0.0, 1.0),
              "1.0"
            )
          )
        ),
        m_MinSizeFlag
        (
          new command::FlagStatic
          (
            "min_size",
            "minimum size of the fragments to keep it in the fragments list",
            command::Parameter
            (
              "min_size",
              "minimum size of the fragments to keep it in the fragments list",
              command::ParameterCheckRanged< size_t>(),
              "2"
            )
          )
        ),
        m_IgnoreScaffoldsWithOpenRingsFlag
        (
          new command::FlagStatic
          (
            "ignore_scaffolds_with_open_rings",
            "whether scaffolds that have broken rings should be ignored",
            command::Parameter
            (
              "yes-no",
              "whether scaffolds that have broken rings should be ignored",
              command::ParameterCheckRanged< bool>(),
              "1"
            )
          )
        ),
        m_IgnoreScaffoldsWithIncompleteRingSystemsFlag
        (
          new command::FlagStatic
          (
            "ignore_scaffolds_with_incomplete_ring_systems",
            "whether scaffolds that have incomplete ring systems should be ignored",
            command::Parameter
            (
              "yes-no",
              "whether scaffolds that have incomplete ring systems should be ignored",
              command::ParameterCheckRanged< bool>(),
              "1"
            )
          )
        ),
        m_BondTypeData
        (
          new command::FlagStatic
          (
            "bond_coloring_scheme",
            "how to color bond data",
            command::Parameter
            (
              "scheme",
              "used to compare bonds",
              command::ParameterCheckSerializable( chemistry::ConfigurationalBondTypeData::DataEnum()),
              "BondOrderInRingOrAromatic"
            )
          )
        ),
        m_OutputFileFlag
        (
          new command::FlagStatic
          (
            "output",
            "filename to output sdf containing scaffolds",
            command::Parameter
            (
              "output_filename",
              "filename to output sdf containing scaffolds"
            )
          )
        )
    {
    }

    //! copy constructor, only copy the flags
    BuildScaffoldLibrary::BuildScaffoldLibrary( const BuildScaffoldLibrary &PARENT) :
          m_EnsembleFileFlag( PARENT.m_EnsembleFileFlag),
          m_SamplingFractionFlag( PARENT.m_SamplingFractionFlag),
          m_MinSizeFlag( PARENT.m_MinSizeFlag),
          m_IgnoreScaffoldsWithOpenRingsFlag( PARENT.m_IgnoreScaffoldsWithOpenRingsFlag),
          m_IgnoreScaffoldsWithIncompleteRingSystemsFlag( PARENT.m_IgnoreScaffoldsWithIncompleteRingSystemsFlag),
          m_BondTypeData( PARENT.m_BondTypeData),
          m_OutputFileFlag( PARENT.m_OutputFileFlag)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BuildScaffoldLibrary
    BuildScaffoldLibrary *BuildScaffoldLibrary::Clone() const
    {
      return new BuildScaffoldLibrary( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &BuildScaffoldLibrary::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> BuildScaffoldLibrary::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // insert all the flags and params
      sp_cmd->AddFlag( m_EnsembleFileFlag);
      sp_cmd->AddFlag( m_SamplingFractionFlag);
      sp_cmd->AddFlag( m_MinSizeFlag);
      sp_cmd->AddFlag( m_IgnoreScaffoldsWithOpenRingsFlag);
      sp_cmd->AddFlag( m_IgnoreScaffoldsWithIncompleteRingSystemsFlag);
      sp_cmd->AddFlag( m_BondTypeData);
      sp_cmd->AddFlag( m_OutputFileFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief determine if a graph, which contains bonds colored by ConfigurationalBondTypeData::e_BondOrderAromaticWithRingness, contains unclosed rings
    //! @param GRAPH the graph, with edges colored by ConfigurationalBondTypeData::e_BondOrderAromaticWithRingness
    //! @param DATA the bond type data data that was used to make the graph
    //! @return true if the graph contains any unclosed rings
    bool BuildScaffoldLibrary::DetermineIfGraphContainsUnclosedRings
    (
      const graph::ConstGraph< size_t, size_t> &GRAPH,
      const chemistry::ConfigurationalBondTypeData::Data &DATA
    )
    {
      // get the minimum value for a ring bond.  The other two ring coloring schemes use 0 to indicate a bond in a ring
      const size_t first_ring_bond_value
      (
        DATA == chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness ? 3 : 0
      );
      for( size_t vertex_id( 0), graph_size( GRAPH.GetSize()); vertex_id < graph_size; ++vertex_id)
      {
        // initialize the list of ring bonds
        storage::List< size_t> ring_bond_ids;
        const storage::Vector< size_t> &neighbors_data( GRAPH.GetNeighborData( vertex_id));
        const storage::Vector< size_t> &neighbors_indices( GRAPH.GetNeighborIndices( vertex_id));
        for
        (
          size_t index_in_neighbors( 0), number_neighbors( neighbors_data.GetSize());
          index_in_neighbors < number_neighbors;
          ++index_in_neighbors
        )
        {
          // bonds are colored 1, 2, or 3 for single, double, triple bonds outside rings, higher number indicate ring membership
          if( neighbors_data( index_in_neighbors) > first_ring_bond_value)
          {
            ring_bond_ids.PushBack( neighbors_indices( index_in_neighbors));
          }
        }
        if( ring_bond_ids.GetSize() == 1)
        {
          return true;
        }
        for
        (
          storage::List< size_t>::const_iterator
          itr_ring_bond_ids( ring_bond_ids.Begin()), itr_ring_bond_ids_end( ring_bond_ids.End());
          itr_ring_bond_ids != itr_ring_bond_ids_end;
          ++itr_ring_bond_ids
        )
        {
          if
          (
            graph::Connectivity::LengthOfSmallestCycleWithEdge( GRAPH, vertex_id, *itr_ring_bond_ids)
            > GRAPH.GetSize()
          )
          {
            return true;
          }
        }
      }
      return false;
    }

    //! @brief determine if a graph, which contains bonds colored by ConfigurationalBondTypeData::e_BondOrderAromaticWithRingness, contains unclosed rings
    //! @param GRAPH the graph, with edges colored by ConfigurationalBondTypeData::e_BondOrderAromaticWithRingness
    //! @param DATA the bond type data data that was used to make the graph
    //! @return true if the graph contains any unclosed rings
    bool BuildScaffoldLibrary::DetermineIfGraphContainsIncompleteRingSystems
    (
      const graph::ConstGraph< size_t, size_t> &GRAPH,
      const chemistry::ConfigurationalBondTypeData::Data &DATA
    )
    {
      // get the minimum value for a ring bond.  The other two ring coloring schemes use 0 to indicate a bond in a ring
      const size_t first_ring_bond_value
      (
        DATA == chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness ? 3 : 0
      );
      for( size_t vertex_id( 0), graph_size( GRAPH.GetSize()); vertex_id < graph_size; ++vertex_id)
      {
        // initialize the list of ring bonds
        storage::List< size_t> ring_bond_ids;
        const storage::Vector< size_t> &neighbors_data( GRAPH.GetNeighborData( vertex_id));
        const storage::Vector< size_t> &neighbors_indices( GRAPH.GetNeighborIndices( vertex_id));
        for
        (
          size_t index_in_neighbors( 0), number_neighbors( neighbors_data.GetSize());
          index_in_neighbors < number_neighbors;
          ++index_in_neighbors
        )
        {
          // bonds are colored 1, 2, or 3 for single, double, triple bonds outside rings, higher number indicate ring membership
          if( neighbors_data( index_in_neighbors) > first_ring_bond_value)
          {
            ring_bond_ids.PushBack( neighbors_indices( index_in_neighbors));
          }
        }
        if( ring_bond_ids.GetSize() == 1)
        {
          return true;
        }
      }
      return false;
    }

    //! @brief make a hash string from a map of size-t's to size-t's
    //! @param MAP the map to make a hashable string from
    //! @return a string containing key,value,key,value pairs
    std::string BuildScaffoldLibrary::MakeHashStringFromMap( const storage::Map< size_t, size_t> &MAP)
    {
      if( MAP.GetSize() == 0)
      {
        return std::string();
      }

      // separator will go between each size_t printed in the map
      // since these are size_ts, it is fine to just use a space
      const char separator( ' ');

      // make a hash beginning with the first element
      std::string hash( util::Format()( MAP.Begin()->first) + separator + util::Format()( MAP.Begin()->second));

      // add all the rest of the elements to the map
      for
      (
        storage::Map< size_t, size_t>::const_iterator itr( ++MAP.Begin()), itr_end( MAP.End());
        itr != itr_end;
        ++itr
      )
      {
        hash += separator + util::Format()( itr->first) + separator + util::Format()( itr->second);
      }

      return hash;
    }

    //! @brief make a hash string from a graph
    //! @param GRAPH a graph
    //! This string can be used as a key to a map that holds graphs with identical hash strings, thus narrowing the
    //! number of graphs that must be searched to determine whether a new scaffold is unique
    //! @return a string that specifies something about the graph that is vertex-invariant
    //!         e.g. does not depend on the ordering of the vertices
    std::string BuildScaffoldLibrary::MakeHashStringFromGraph( const graph::ConstGraph< size_t, size_t> &GRAPH)
    {
      // make maps containing counts of each vertex and edge type
      const storage::Map< size_t, size_t> &vertex_counts( GRAPH.GetVertexTypeCountMap());
      const storage::Map< size_t, size_t> &edge_counts( GRAPH.GetEdgeTypeCountMap());

      // make hash strings from each map, separated by "$"
      return MakeHashStringFromMap( vertex_counts) + "$" + MakeHashStringFromMap( edge_counts);
    }

    //! @brief get the next pair to compare
    storage::Pair< size_t, size_t> BuildScaffoldLibrary::GetNextPairToCompare() const
    {
      m_GetNextPairMutex.Lock();

      // get the size of the ensemble
      static const size_t ensemble_size( m_EnsembleSimpleGraphs.GetSize());

      // get the fraction of pairs to sample
      static const double fraction( m_SamplingFractionFlag->GetFirstParameter()->GetNumericalValue< double>());

      static math::Range< double> random_range( 0.0, 1.0);

      while( m_LastAssignedMoleculeIds.First() < ensemble_size)
      {
        while( m_LastAssignedMoleculeIds.Second() < m_LastAssignedMoleculeIds.First())
        {
          if( random::GetGlobalRandom().Double( random_range) <= fraction)
          {
            storage::Pair< size_t, size_t> next_pair( m_LastAssignedMoleculeIds);
            ++m_LastAssignedMoleculeIds.Second();
            m_GetNextPairMutex.Unlock();
            return next_pair;
          }
          ++m_LastAssignedMoleculeIds.Second();
        }
        ++m_LastAssignedMoleculeIds.First();
        m_LastAssignedMoleculeIds.Second() = 0;
      }

      // return a pair that is out of bounds, at which point the thread will stop
      storage::Pair< size_t, size_t> next_pair( m_LastAssignedMoleculeIds);
      m_GetNextPairMutex.Unlock();
      return next_pair;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int BuildScaffoldLibrary::Main() const
    {
      // read in ensemble
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_EnsembleFileFlag->GetFirstParameter()->GetValue());
      // remove any hydrogens because they slow down the scaffold search and are unnecessary
      chemistry::FragmentEnsemble ensemble( input, sdf::e_Remove);
      io::File::CloseClearFStream( input);

      // store the size of the ensemble
      const size_t ensemble_size( ensemble.GetSize());

      if( ensemble_size <= 1)
      {
        BCL_MessageCrt( "Cannot find common scaffolds without at least two molecules");
        return 1;
      }

      // Get the coloring scheme to be used for bonds
      m_BondColoringScheme =
        chemistry::ConfigurationalBondTypeData::DataEnum( m_BondTypeData->GetFirstParameter()->GetValue());

      // determine whether we should skip graphs that have open rings
      const bool skip_graphs_with_open_rings( m_IgnoreScaffoldsWithOpenRingsFlag->GetFirstParameter()->GetNumericalValue< bool>());
      const bool skip_graphs_with_incomplete_rings( m_IgnoreScaffoldsWithIncompleteRingSystemsFlag->GetFirstParameter()->GetNumericalValue< bool>());

      if
      (
        ( skip_graphs_with_open_rings || skip_graphs_with_incomplete_rings)
        && m_BondColoringScheme != chemistry::ConfigurationalBondTypeData::e_IsInRing
        && m_BondColoringScheme != chemistry::ConfigurationalBondTypeData::e_BondOrderInRingOrAromatic
        && m_BondColoringScheme != chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
      )
      {
        BCL_MessageCrt
        (
          "Cannot skip bonds with open rings unless bond comparison type is one of\n"
          "IsInRing, BondOrderInRingOrAromatic, or BondOrderOrAromaticWithRingness"
        );
      }

      // For each small molecule, instantiate a graphical representation based on atom type and bond order/aromatic
      // store the graphs in vectors
      m_EnsembleSimpleGraphs.Reset();
      m_EnsembleSimpleGraphs.AllocateMemory( ensemble_size);
      m_EnsembleGraphs.Reset();
      m_EnsembleGraphs.AllocateMemory( ensemble_size);
      chemistry::ConformationGraphConverter graph_converter
      (
        chemistry::ConformationGraphConverter::e_AtomType,
        m_BondColoringScheme
      );
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        m_EnsembleSimpleGraphs.PushBack( graph_converter( *itr));
        m_EnsembleGraphs.PushBack( chemistry::ConformationGraphConverter::CreateGraphWithAtoms( *itr));
      }

      // reset # of pairs seen
      m_NumberOfScaffolds = 0;
      m_NumberPairsToConsider = 0;
      m_LastAssignedMoleculeIds = storage::Pair< size_t, size_t>( 0, 0);

      const size_t initial_seed( random::GetGlobalRandom().GetSeed());

      // get the fraction of pairs to sample
      const double fraction( m_SamplingFractionFlag->GetFirstParameter()->GetNumericalValue< double>());

      // allocate pairs equally among all processor's vectors
      for( size_t index_a( 0); index_a < ensemble_size; ++index_a)
      {
        for( size_t index_b( 0); index_b < index_a; ++index_b)
        {
          if( random::GetGlobalRandom().Double( math::Range< double>( 0.0, 1.0)) <= fraction)
          {
            ++m_NumberPairsToConsider;
          }
        }
      }

      random::GetGlobalRandom().SetSeed( initial_seed);

      // open the output file so that we can write out scaffolds as we go
      io::File::MustOpenOFStream( m_Output, m_OutputFileFlag->GetFirstParameter()->GetValue());

      // get the number of processors
      const size_t n_processors( sched::GetNumberCPUs());

      // create a vector to hold the jobs
      util::ShPtrVector< sched::JobInterface> jobs;

      // make the job vector big enough to hold all the jobs
      jobs.AllocateMemory( n_processors);
      const size_t group_id( 0);
      for( size_t processor_number( 0); processor_number < n_processors; ++processor_number)
      {
        // create the job
        jobs.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::ThunkJob< BuildScaffoldLibrary, void>
            (
              group_id,
              *this,
              &BuildScaffoldLibrary::RunThread,
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );

        // submit it to the scheduler
        sched::GetScheduler().SubmitJob( jobs.LastElement());
      }

      // join all the jobs
      for( size_t processor_number( 0); processor_number < n_processors; ++processor_number)
      {
        sched::GetScheduler().Join( jobs( processor_number));
      }

      BCL_MessageStd
      (
        "Saved the " + util::Format()( m_NumberOfScaffolds)
        + " scaffolds found to " + m_OutputFileFlag->GetFirstParameter()->GetValue()
      );

      // close the file stream
      io::File::CloseClearFStream( m_Output);

      return 0;
    }

    //! @brief test whether an isomorphism is large enough
    //! @param SUBGRAPH_ISOMORPHISM_SIZE the size of the isomorphism
    //! @return true if the subgraph isomorphism size is larger than or equal to the m_MinSizeFlag
    bool BuildScaffoldLibrary::IsIsomorphismLargerThanMinSize( const size_t &SUBGRAPH_ISOMORPHISM_SIZE) const
    {
      static const size_t s_MinIsomorphismSize( m_MinSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>());
      return SUBGRAPH_ISOMORPHISM_SIZE >= s_MinIsomorphismSize;
    }

    //! @brief increment the count of pairs we have checked and notify the user if we are far enough along
    void BuildScaffoldLibrary::IncrementPairsCompared() const
    {
      // lock the writing mutex so that other threads cannot increment pairs or m_NumberOfScaffolds
      m_WritingMutex.Lock();
      static size_t pairs( 0);

      // notify every 500 pairs or every time that the estimated pairs to compare would show an increase
      static const size_t s_PairsPerNotification( std::max( size_t( 500), m_NumberPairsToConsider / 10000));

      // notify the user whenever pairs is a multiple of s_PairsPerNotification
      if( !( ++pairs % s_PairsPerNotification))
      {
        BCL_MessageStd
        (
          "% complete: "
          + util::Format().FFP( 2)( 100.0 * double( pairs) / double( m_NumberPairsToConsider))
          + " # scaffolds: " + util::Format()( m_NumberOfScaffolds)
        );
      }

      // unlock the mutex so that other threads can write out small molecules or increment pairs
      m_WritingMutex.Unlock();
    }

    //! @brief compare the molecules given by the indices in a vector
    void BuildScaffoldLibrary::RunThread() const
    {
      // this object will hold the isomorphism between any graphs that are compared (vertices colored by atom type index)
      graph::CommonSubgraphIsomorphism< size_t, size_t>
      isomorphism( graph::CommonSubgraphIsomorphismBase::e_Connected);

      // this object will hold the isomorphism between the graphs (vertices colored by si-ptr to atom)
      graph::CommonSubgraphIsomorphism< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t>
      isomorphism_with_atoms( graph::CommonSubgraphIsomorphismBase::e_Connected);

      // store the minimum size for an isomorphism to be interesting
      const size_t min_size( m_MinSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>());

      // determine whether we should skip graphs that have open rings
      const bool skip_graphs_with_open_rings
      (
        m_IgnoreScaffoldsWithOpenRingsFlag->GetFirstParameter()->GetNumericalValue< bool>()
        &&
        (
          m_BondColoringScheme == chemistry::ConfigurationalBondTypeData::e_IsInRing
          || m_BondColoringScheme == chemistry::ConfigurationalBondTypeData::e_BondOrderInRingOrAromatic
          || m_BondColoringScheme == chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
        )
      );
      // determine whether we should skip graphs that have incomplete rings
      const bool skip_graphs_with_incomplete_rings
      (
        m_IgnoreScaffoldsWithIncompleteRingSystemsFlag->GetFirstParameter()->GetNumericalValue< bool>()
        &&
        (
          m_BondColoringScheme == chemistry::ConfigurationalBondTypeData::e_IsInRing
          || m_BondColoringScheme == chemistry::ConfigurationalBondTypeData::e_BondOrderInRingOrAromatic
          || m_BondColoringScheme == chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
        )
      );

      // get the name of the file that all the molecules are from
      const std::string ensemble_filename( m_EnsembleFileFlag->GetFirstParameter()->GetValue());

      // store the size of the ensemble
      const size_t ensemble_size( m_EnsembleGraphs.GetSize());

      // keep track of what molecule # we are on
      for
      (
        storage::Pair< size_t, size_t> pair_to_compare( GetNextPairToCompare());
        pair_to_compare.First() < ensemble_size;
        pair_to_compare = GetNextPairToCompare()
      )
      {
        // get the indices of the molecules out of the pair
        const size_t mols_a_indice( pair_to_compare.First()), mols_b_indice( pair_to_compare.Second());

        // increment the # of pairs compared so far
        IncrementPairsCompared();

        // set up the isomorphisms with pointers to the graphs
        // copying the graphs is fairly expensive, so we avoid it whenever possible
        isomorphism.SetGraphA
        (
          util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &m_EnsembleSimpleGraphs( mols_a_indice), false)
        );
        isomorphism.SetGraphB
        (
          util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &m_EnsembleSimpleGraphs( mols_b_indice), false)
        );

        // get the actual isomorphism, using estimated upper bounds on its size
        isomorphism.FindIsomorphism( isomorphism.EstimateUpperBounds(), min_size);

        // if the isomorphism is too small, ignore it
        if( !IsIsomorphismLargerThanMinSize( isomorphism.GetIsomorphism().GetSize()))
        {
          continue;
        }

        // get the actual subgraph from graph a (this choice is arbitrary)
        graph::ConstGraph< size_t, size_t> subgraph( isomorphism.GetSubgraphIsomorphismsOfGraphA().FirstElement().ToGraph());

        // skip graphs with open rings, if the appropriate flag was set
        if( skip_graphs_with_open_rings && DetermineIfGraphContainsUnclosedRings( subgraph, m_BondColoringScheme))
        {
          continue;
        }
        else if( skip_graphs_with_incomplete_rings && DetermineIfGraphContainsIncompleteRingSystems( subgraph, m_BondColoringScheme))
        {
          continue;
        }

        // set isomorphism_with_atoms to have the correct graphs and isomorphism
        isomorphism_with_atoms.SetGraphA
        (
          util::OwnPtr< graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> >
          (
            &m_EnsembleGraphs( mols_a_indice),
            false
          )
        );
        isomorphism_with_atoms.SetGraphB
        (
          util::OwnPtr< graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> >
          (
            &m_EnsembleGraphs( mols_b_indice),
            false
          )
        );
        isomorphism_with_atoms.SetIsomorphisms( isomorphism.GetIsomorphisms());

        // set the isomorphism up to search for the subgraph in similar scaffolds, give it a non-owning pointer
        // to avoid costly copies
        isomorphism.SetGraphA( util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &subgraph, false));

        // set this to false if we find that the scaffold is non-unique
        bool scaffold_is_unique( true);

        // get a hash-string from the graph to use as a key in the unique scaffolds map
        // this way, we only need to check subgraph against other graphs with the same hash (which should be very few)
        const std::string subgraph_isomorphism_hash( MakeHashStringFromGraph( subgraph));

        // so that we can read/write to the unique scaffolds hash map, we need to briefly lock the mutex
        m_UniqueScaffoldsMutex.Lock();

        // get an iterator to the position in the map with the same hash
        storage::List< graph::ConstGraph< size_t, size_t> > &similar_scaffolds( m_UniqueScaffolds[ subgraph_isomorphism_hash]);
        sched::Mutex &list_mutex( m_ScaffoldListMuteces[ subgraph_isomorphism_hash]);

        // unlock the mutex, since we can now lock a mutex for just the scaffolds with the same hash
        m_UniqueScaffoldsMutex.Unlock();

        // lock the mutex for the list of scaffolds with the same hash
        list_mutex.Lock();

        // examine scaffolds with the same hash in the map to test whether this scaffold is unique
        for
        (
          storage::List< graph::ConstGraph< size_t, size_t> >::iterator
          itr_similar_scaffolds( similar_scaffolds.Begin()),
          itr_similar_scaffolds_end( similar_scaffolds.End());
          itr_similar_scaffolds != itr_similar_scaffolds_end;
          ++itr_similar_scaffolds
        )
        {
          // get the matching vertices in the similar scaffold to the current subgraph
          const storage::Vector< storage::Vector< size_t> > matching_vertices
          (
            graph::CSISubstructure::GetVertexMatchingMatrixForSubgraph( subgraph, *itr_similar_scaffolds)
          );

          // GetVertexMatchingMatrixForSubgraph stops as soon as any vertices do not have a compatible vertex in
          // the larger graph, so if the last element is empty, then the graphs cannot be equal
          if( matching_vertices.LastElement().IsEmpty())
          {
            continue;
          }

          // it is very likely at this point that the scaffold and the subgraph are the same
          // Nevertheless, we need to check with the isomorphism to be sure

          // make a non-owning pointer to the existing scaffold
          util::OwnPtr< graph::ConstGraph< size_t, size_t> > op_known_scaffold( &*itr_similar_scaffolds, false);

          // set up isomorphism to check whether the known scaffold is identical to subgraph
          isomorphism.SetGraphB( op_known_scaffold);

          // perform the isomorphism search
          isomorphism.FindIsomorphism( subgraph.GetSize(), subgraph.GetSize(), matching_vertices);

          // if a complete isomorphism was found, the scaffold is not unique
          if( isomorphism.GetIsomorphism().GetSize() == subgraph.GetSize())
          {
            scaffold_is_unique = false;

            // break out of this loop since we now know that the scaffold was seen before
            break;
          }
        }

        // if the scaffold was seen before, ignore it
        if( !scaffold_is_unique)
        {
          // unlock the mutex for the molecules with the same hash, since we're going on to the next pair of molecules
          list_mutex.Unlock();
          continue;
        }

        // a novel scaffold was found

        // add it to the list of similar scaffolds
        similar_scaffolds.PushBack( subgraph);

        // unlock the mutex for the molecules with the same hash
        list_mutex.Unlock();

        // for a name for the new scaffold using the molecule #'s from the ensemble file
        const std::string new_name
        (
          "Scaffold produced from molecule #s "
          + util::Format()( mols_a_indice) + " & " + util::Format()( mols_b_indice)
          + "\nfrom file\n" + ensemble_filename
        );

        // convert the scaffold graph back into a small molecule and add properties of interest
        util::SiPtr< chemistry::FragmentComplete> sp_new_scaffold
        (
          new chemistry::FragmentComplete
          (
            chemistry::ConformationGraphConverter::CreateAtomsFromGraph
            (
              isomorphism_with_atoms.GetSubgraphIsomorphismsOfGraphA().FirstElement().ToGraph()
            ),
            new_name
          )
        );

        //sp_new_scaffold->GetProperties().SetMDLProperty( "HasOpenBonds", util::Format()( graph_had_open_rings));

        // lock the writing mutex so that we can increment m_NumberOfScaffolds and write the scaffold out to a file
        m_WritingMutex.Lock();

        // set the bcl scaffold # to the index of the molecule in scaffolds
        sp_new_scaffold->StoreProperty( "BclScaffoldNumber", util::Format()( m_NumberOfScaffolds));

        // increment the number of scaffolds seen
        ++m_NumberOfScaffolds;

        // write the new scaffold out to a file
        sp_new_scaffold->WriteMDL( m_Output);

        // unlock the writing mutex, now that we are done writing to the file
        m_WritingMutex.Unlock();
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BuildScaffoldLibrary::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &BuildScaffoldLibrary::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
