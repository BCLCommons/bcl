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

#ifndef BCL_ASSEMBLE_AA_NEIGHBOR_LIST_CONTAINER_GENERATOR_SSE_H_
#define BCL_ASSEMBLE_AA_NEIGHBOR_LIST_CONTAINER_GENERATOR_SSE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AANeighborListContainerGeneratorSSE
    //! @brief BinaryFunction derived class, that generates the AANeighborlistContainer from two secondary structure elements
    //! @details
    //!
    //! @see @link example_assemble_aa_neighbor_list_container_cenerator_sse_pair.cpp @endlink
    //! @author woetzen
    //! @date Apr 21, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AANeighborListContainerGeneratorSSE :
      public math::FunctionInterfaceSerializable< SSE, AANeighborListContainer>
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

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief construct from distance cutoff and minimal sequence separation
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      AANeighborListContainerGeneratorSSE
      (
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN
      );

    public:

      //! @brief Clone function
      //! @return pointer to new AANeighborListContainerGeneratorSSE
      AANeighborListContainerGeneratorSSE *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief AANeighborlist from two secondary structure elements
      //! @param ELEMENT secondary structure element
      //! @return the AANeighborlistContainer constructed from all amino in the SSE
      AANeighborListContainer operator()( const SSE &ELEMENT) const;

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

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return a shptr to NeighborListGenerator class
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      //! @param CACHED wrap the generator in a cache object
      //! @return ShPtr to generator
      static util::ShPtr< math::FunctionInterfaceSerializable< SSE, AANeighborListContainer> > AANeighborListGenerator
      (
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN,
        const bool CACHED
      );

    }; // class AANeighborListContainerGeneratorSSE

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_AA_NEIGHBOR_LIST_CONTAINER_GENERATOR_SSE_H_ 
