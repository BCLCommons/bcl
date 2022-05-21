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
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_sse.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_sse_pair.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_pruner.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_function_adapter.h"
#include "math/bcl_math_function_cached.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AANeighborListContainerGeneratorProteinModel::s_Instance
    (
      GetObjectInstances().AddInstance( new AANeighborListContainerGeneratorProteinModel( 0.0, 0, true, false))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from SSE pair neighbor list container generator
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    AANeighborListContainerGeneratorProteinModel::AANeighborListContainerGeneratorProteinModel
    (
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN,
      const bool CACHED
    ) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_MinimalSequenceSeparation( MIN_SEQ_SEPARATION),
      m_ConsiderDifferentChain( CONSIDER_DIFFERENT_CHAIN),
      m_Cached( CACHED)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AANeighborListContainerGeneratorProteinModel
    AANeighborListContainerGeneratorProteinModel *AANeighborListContainerGeneratorProteinModel::Clone() const
    {
      return new AANeighborListContainerGeneratorProteinModel( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AANeighborListContainerGeneratorProteinModel::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief generate AANeighborListContainer for protein model
    //! @param PROTEIN_MODEL the protein model
    //! @return AANeighborListContainer for ProteinModel
    AANeighborListContainer AANeighborListContainerGeneratorProteinModel::operator()( const ProteinModel &PROTEIN_MODEL) const
    {
      return
        AANeighborListContainer
        (
          PROTEIN_MODEL.GetAminoAcids(), m_DistanceCutoff, m_MinimalSequenceSeparation, m_ConsiderDifferentChain
        );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AANeighborListContainerGeneratorProteinModel::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceCutoff               , ISTREAM);
      io::Serialize::Read( m_MinimalSequenceSeparation    , ISTREAM);
      io::Serialize::Read( m_ConsiderDifferentChain       , ISTREAM);
      io::Serialize::Read( m_Cached,                        ISTREAM);
      io::Serialize::Read( m_AANeigborListGeneratorSSEPair, ISTREAM);
      io::Serialize::Read( m_AANeigborListGeneratorSSE    , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AANeighborListContainerGeneratorProteinModel::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceCutoff               , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_MinimalSequenceSeparation    , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ConsiderDifferentChain       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Cached                       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AANeigborListGeneratorSSEPair, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AANeigborListGeneratorSSE    , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return a shptr to NeighborListGenerator class
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    //! @param CACHED wrap the generator in a cache object
    //! @return ShPtr to generator
    util::ShPtr< math::FunctionInterfaceSerializable< ProteinModel, AANeighborListContainer> >
    AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
    (
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN,
      const bool CACHED
    )
    {
      static double s_max_distance_cutoff_for_caching( 8);
      static size_t s_min_seq_separation_for_caching( 0);
      // no caching desired or distance is so large as to make caching futile due to the amount of pruning necessary for
      // common applications
      if
      (
        !CACHED
        || DISTANCE_CUTOFF > s_max_distance_cutoff_for_caching
        || MIN_SEQ_SEPARATION < s_min_seq_separation_for_caching
      )
      {
        // just create the generator
        return util::ShPtr< math::FunctionInterfaceSerializable< ProteinModel, AANeighborListContainer> >
        (
          new AANeighborListContainerGeneratorProteinModel
          (
            DISTANCE_CUTOFF, MIN_SEQ_SEPARATION, CONSIDER_DIFFERENT_CHAIN, CACHED
          )
        );
      }

      // caching desired, create a generator that make the maximal possible neighbor list
      typedef storage::Map< std::string, util::ShPtr< math::FunctionCached< ProteinModel, AANeighborListContainer> > > PrunersMapType;

      // initialize static map to store the pruners
      static PrunersMapType s_pruners;

      // construct the identifier for the requested pruner
      const std::string identifier
      (
        util::Format()( DISTANCE_CUTOFF) + '_' +
        util::Format()( MIN_SEQ_SEPARATION) + '_' +
        util::Format()( CONSIDER_DIFFERENT_CHAIN)
      );

      // search for this identifier in the map to see if such a pruner was previously constructed
      const PrunersMapType::const_iterator pruner_itr( s_pruners.Find( identifier));

      // if found, then return this pruner
      if( pruner_itr != s_pruners.End())
      {
        return pruner_itr->second;
      }

      // caching desired, create a static generator that make the maximal possible neighbor list
      static util::ShPtr< math::FunctionCached< ProteinModel, AANeighborListContainer> >
        s_complete_neighbor_list_cached;
      // if not initialize yet, then initialize it
      if( !s_complete_neighbor_list_cached.IsDefined())
      {
        // construct generator
        const util::ShPtr< math::FunctionInterfaceSerializable< ProteinModel, AANeighborListContainer> >
          complete_neighbor_list
          (
            new AANeighborListContainerGeneratorProteinModel
            (
              s_max_distance_cutoff_for_caching,
              s_min_seq_separation_for_caching,
              true,
              true
            )
          );
        // wrap it in cache
        s_complete_neighbor_list_cached = util::ShPtr< math::FunctionCached< ProteinModel, AANeighborListContainer> >
        (
          new math::FunctionCached< ProteinModel, AANeighborListContainer>
          (
            complete_neighbor_list,
            &ProteinModel::GetDestructorSignal,
            &ProteinModel::GetChangeSignal
          )
        );
      }

      if
      (
        DISTANCE_CUTOFF == s_max_distance_cutoff_for_caching
        && MIN_SEQ_SEPARATION == s_min_seq_separation_for_caching
        && CONSIDER_DIFFERENT_CHAIN
      )
      {
        s_pruners[ identifier] = s_complete_neighbor_list_cached;
        return s_complete_neighbor_list_cached;
      }

      // otherwise construct the pruner
      const util::ShPtr< math::FunctionInterfaceSerializable< AANeighborListContainer, AANeighborListContainer> > neighborlist_pruner
      (
        new AANeighborListContainerPruner
        (
          DISTANCE_CUTOFF, MIN_SEQ_SEPARATION, CONSIDER_DIFFERENT_CHAIN
        )
      );
      // wrap cache around it
      util::ShPtr< math::FunctionCached< ProteinModel, AANeighborListContainer> > neighborlist_pruner_cached
      (
        new math::FunctionCached< ProteinModel, AANeighborListContainer>
        (
          util::ShPtr< math::FunctionInterfaceSerializable< ProteinModel, AANeighborListContainer> >
          (
            new math::FunctionAdapter< ProteinModel, AANeighborListContainer, AANeighborListContainer>
            (
              s_complete_neighbor_list_cached, neighborlist_pruner
            )
          ),
          &ProteinModel::GetDestructorSignal,
          &ProteinModel::GetChangeSignal
        )
      );
      // insert it into the map and return it
      s_pruners[ identifier] = neighborlist_pruner_cached;
      return neighborlist_pruner_cached;
    }

  } // namespace assemble
} // namespace bcl
