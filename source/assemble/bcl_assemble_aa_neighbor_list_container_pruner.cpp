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
#include "assemble/bcl_assemble_aa_neighbor_list_container_pruner.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AANeighborListContainerPruner::s_Instance
    (
      GetObjectInstances().AddInstance( new AANeighborListContainerPruner( 0.0, 0, true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from distance cutoff and minimal sequence separation
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    AANeighborListContainerPruner::AANeighborListContainerPruner
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
    //! @return pointer to new AANeighborListContainerPruner
    AANeighborListContainerPruner *AANeighborListContainerPruner::Clone() const
    {
      return new AANeighborListContainerPruner( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AANeighborListContainerPruner::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that prunes a given neighbor list container to a reduced one
    //! @param CONTAINER the container to be pruned
    //! @return pruned AANeighborListContainer
    AANeighborListContainer AANeighborListContainerPruner::operator()( const AANeighborListContainer &CONTAINER) const
    {
      return AANeighborListContainer( CONTAINER, m_DistanceCutoff, m_MinimalSequenceSeparation, m_ConsiderDifferentChain);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AANeighborListContainerPruner::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AANeighborListContainerPruner::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
