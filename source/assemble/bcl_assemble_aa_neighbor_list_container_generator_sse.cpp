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
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_sse.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_pruner.h"
#include "assemble/bcl_assemble_sse.h"
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
    const util::SiPtr< const util::ObjectInterface> AANeighborListContainerGeneratorSSE::s_Instance
    (
      GetObjectInstances().AddInstance( new AANeighborListContainerGeneratorSSE( 0.0, 0, true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from distance cutoff and minimal sequence separation
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    AANeighborListContainerGeneratorSSE::AANeighborListContainerGeneratorSSE
    (
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN
    ) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_MinimalSequenceSeparation( MIN_SEQ_SEPARATION),
      m_ConsiderDifferentChain( CONSIDER_DIFFERENT_CHAIN)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AANeighborListContainerGeneratorSSE
    AANeighborListContainerGeneratorSSE *AANeighborListContainerGeneratorSSE::Clone() const
    {
      return new AANeighborListContainerGeneratorSSE( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AANeighborListContainerGeneratorSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief AANeighborlist from two secondary structure elements
    //! @param ELEMENT secondary structure element
    //! @return the AANeighborlistContainer constructed from all amino in the SSE
    AANeighborListContainer AANeighborListContainerGeneratorSSE::operator()( const SSE &ELEMENT) const
    {
      // calculate neighbor list container and return
      return
        AANeighborListContainer
        (
          ELEMENT.GetMembers(), m_DistanceCutoff, m_MinimalSequenceSeparation, m_ConsiderDifferentChain
        );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AANeighborListContainerGeneratorSSE::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceCutoff           , ISTREAM);
      io::Serialize::Read( m_MinimalSequenceSeparation, ISTREAM);
      io::Serialize::Read( m_ConsiderDifferentChain   , ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AANeighborListContainerGeneratorSSE::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceCutoff           , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_MinimalSequenceSeparation, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ConsiderDifferentChain   , OSTREAM, INDENT);

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
    util::ShPtr< math::FunctionInterfaceSerializable< SSE, AANeighborListContainer> > AANeighborListContainerGeneratorSSE::AANeighborListGenerator
    (
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN,
      const bool CACHED
    )
    {
      static size_t s_max_distance_cutoff_for_caching( 12.0);
      // no caching desired or distance is so large as to make caching futile due to the amount of pruning necessary for
      // common applications
      if( !CACHED || DISTANCE_CUTOFF > s_max_distance_cutoff_for_caching)
      {
        // just create the generator
        return util::ShPtr< math::FunctionInterfaceSerializable< SSE, AANeighborListContainer> >
        (
          new AANeighborListContainerGeneratorSSE( DISTANCE_CUTOFF, MIN_SEQ_SEPARATION, CONSIDER_DIFFERENT_CHAIN)
        );
      }

      // caching desired, create a generator that make the maximal possible neighbor list
      static util::ShPtr< math::FunctionCached< SSE, AANeighborListContainer> > s_complete_neighbor_list_cached;
      if( !s_complete_neighbor_list_cached.IsDefined())
      {
        const util::ShPtr< math::FunctionInterfaceSerializable< SSE, AANeighborListContainer> > complete_neighbor_list
        (
          new AANeighborListContainerGeneratorSSE( s_max_distance_cutoff_for_caching, 0, true)
        );
        s_complete_neighbor_list_cached = util::ShPtr< math::FunctionCached< SSE, AANeighborListContainer> >
        (
          new math::FunctionCached< SSE, AANeighborListContainer>( complete_neighbor_list, &SSE::GetDestructorSignal)
        );
        s_complete_neighbor_list_cached->AddSignalHandlerForArgument( &SSE::GetCoordinateChangeSignal);
      }

      // use a pruner, to narrow the neighbor list
      typedef storage::Map< std::string, util::ShPtr< math::FunctionCached< SSE, AANeighborListContainer> > > PrunersMapType;
      static PrunersMapType s_pruners;
      const std::string identifier
      (
        util::Format()( DISTANCE_CUTOFF) + '_' +
        util::Format()( MIN_SEQ_SEPARATION) + '_' +
        util::Format()( CONSIDER_DIFFERENT_CHAIN)
      );
      const PrunersMapType::const_iterator pruner_itr( s_pruners.Find( identifier));
      if( pruner_itr != s_pruners.End())
      {
        return pruner_itr->second;
      }

      const util::ShPtr< math::FunctionInterfaceSerializable< SSE, AANeighborListContainer> > neighborlist_pruner
      (
        new math::FunctionAdapter< SSE, AANeighborListContainer, AANeighborListContainer>
        (
          s_complete_neighbor_list_cached,
          util::ShPtr< math::FunctionInterfaceSerializable< AANeighborListContainer, AANeighborListContainer> >
          (
            new AANeighborListContainerPruner( DISTANCE_CUTOFF, MIN_SEQ_SEPARATION, CONSIDER_DIFFERENT_CHAIN)
          )
        )
      );

      util::ShPtr< math::FunctionCached< SSE, AANeighborListContainer> > neighborlist_pruner_cached
      (
        new math::FunctionCached< SSE, AANeighborListContainer>( neighborlist_pruner, &SSE::GetDestructorSignal)
      );
      neighborlist_pruner_cached->AddSignalHandlerForArgument( &SSE::GetCoordinateChangeSignal);
      s_pruners[ identifier] = neighborlist_pruner_cached;

      return neighborlist_pruner_cached;
    }

  } // namespace assemble
} // namespace bcl
