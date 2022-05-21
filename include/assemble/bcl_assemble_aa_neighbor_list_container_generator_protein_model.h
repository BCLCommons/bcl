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

#ifndef BCL_ASSEMBLE_AA_NEIGHBOR_LIST_CONTAINER_GENERATOR_PROTEIN_MODEL_H_
#define BCL_ASSEMBLE_AA_NEIGHBOR_LIST_CONTAINER_GENERATOR_PROTEIN_MODEL_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AANeighborListContainerGeneratorProteinModel
    //! @brief calculate the neighbor list container for a ProteinModel
    //! @details
    //!
    //! @see @link example_assemble_aa_neighbor_list_container_generator_protein_model.cpp @endlink
    //! @author woetzen
    //! @date Apr 21, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AANeighborListContainerGeneratorProteinModel :
      public math::FunctionInterfaceSerializable< ProteinModel, AANeighborListContainer>
    {

    private:

    //////////
    // data //
    //////////

      //! distance cutoff for neighbor list container
      double m_DistanceCutoff;

      //! minimal sequence separation for neighbor list container
      size_t m_MinimalSequenceSeparation;

      //! consider different chain ids for neighbor list container
      bool m_ConsiderDifferentChain;

      //! boolean whether results are cached or not
      bool m_Cached;

      //! @brief function that calculates the neighbor list for a pair of sses
      util::ShPtr< math::BinaryFunctionInterface< SSE, SSE, AANeighborListContainer> > m_AANeigborListGeneratorSSEPair;

      //! @brief function that calculates the neighbor list for a sse
      util::ShPtr< math::FunctionInterfaceSerializable< SSE, AANeighborListContainer> > m_AANeigborListGeneratorSSE;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from SSE pair neighbor list container generator
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      //! @param CACHED
      AANeighborListContainerGeneratorProteinModel
      (
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN,
        const bool CACHED
      );

    public:

      //! @brief Clone function
      //! @return pointer to new AANeighborListContainerGeneratorProteinModel
      AANeighborListContainerGeneratorProteinModel *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief generate AANeighborListContainer for protein model
      //! @param PROTEIN_MODEL the protein model
      //! @return AANeighborListContainer for ProteinModel
      AANeighborListContainer operator()( const ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief return a shptr to NeighborListGenerator class
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      //! @param CACHED wrap the generator in a cache object
      //! @return ShPtr to generator
      static util::ShPtr< math::FunctionInterfaceSerializable< ProteinModel, AANeighborListContainer> >
      AANeighborListGenerator
      (
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN,
        const bool CACHED
      );

    }; // class AANeighborListContainerGeneratorProteinModel

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_AA_NEIGHBOR_LIST_CONTAINER_GENERATOR_PROTEIN_MODEL_H_
