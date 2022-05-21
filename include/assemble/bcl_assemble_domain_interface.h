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

#ifndef BCL_ASSEMBLE_DOMAIN_INTERFACE_H_
#define BCL_ASSEMBLE_DOMAIN_INTERFACE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse.h"
#include "coord/bcl_coord_movable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DomainInterface
    //! @brief Interface from which all Domain classes can derived from including Sheet, Domain, etc.
    //! @details This interface provides a framework to have a collection of SSEs with a general orientation information
    //! The classes are not required to have ownership of the SSEs, they can have just SiPtrs and the SSEs can't be
    //! changed through the interface, this makes sure all the SSE changes has to go through ProteinModel
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Jul 20, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DomainInterface :
      public coord::MovableInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief get SiPtrVector of SSEs
      //! @return SiPtrVector of SSEs
      virtual util::SiPtrVector< const SSE> GetSSEs() const = 0;

      //! @brief returns all SSEs in domain of given SSTYPE in a util::SiPtrVector
      //! @param SS_TYPE specific SSTYPE
      //! @return all SSEs in domain of given SSTYPE in a util::SiPtrVector
      virtual util::SiPtrVector< const SSE> GetSSEs( const biol::SSType &SS_TYPE) const;

      //! @brief returns all SSEs in domain of given SSTYPEs in a util::SiPtrVector
      //! @param SS_TYPES set of sstypes
      //! @return all SSEs in domain of given SSTYPEs in a util::SiPtrVector
      virtual util::SiPtrVector< const SSE> GetSSEs( const storage::Set< biol::SSType> &SS_TYPES) const;

      //! @brief return total number of sses
      //! @return total number of sses
      virtual size_t GetNumberSSEs() const;

      //! @brief returns the number of amino acids in the chain
      //! @return the number of amino acids in the chain
      virtual size_t GetNumberAAs() const = 0;

      //! @brief return numbers of SSE of specified SSTYPE
      //! @param SS_TYPE specific SSTYPE
      //! @return number of SSE of specified SSTYPE
      virtual size_t GetNumberSSEs( const biol::SSType &SS_TYPE) const;

      //! @brief returns all atoms in domain as SiPtrVector
      //! @return all atoms in domain as SiPtrVector
      virtual util::SiPtrVector< const biol::Atom> GetAtoms() const;

      //! @brief returns all atoms of specified ATOM_TYPES in domain as SiPtrVector
      //! @param ATOM_TYPES Set of AtomTypes of interest
      //! @return all atoms of specified types in domain as SiPtrVector
      virtual util::SiPtrVector< const biol::Atom> GetAtoms( const storage::Set< biol::AtomType> &ATOM_TYPES) const;

      //! @brief returns coordinates for all atoms in domain as SiPtrVector
      //! @return coordinates for all atoms in domain as SiPtrVector
      virtual util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates() const;

      //! @brief returns coordinates for all atoms of specified ATOM_TYPES in domain as SiPtrVector
      //! @param ATOM_TYPES Set of AtomTypes of interest
      //! @returns coordinates for all atoms of specified types in domain as SiPtrVector
      virtual util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates
      (
        const storage::Set< biol::AtomType> &ATOM_TYPES
      ) const;

      //! @brief returns the center of the domain
      //! @return the center of the domain
      virtual linal::Vector3D GetCenter() const;

      //! @brief return all the amino acids in all SSEs
      //! @return all the amino acids in all SSEs
      virtual util::SiPtrVector< const biol::AABase> GetAminoAcids() const;

      //! @brief create and return SSEGeometries for all SSEs in this Domain
      //! @return ShPtrVector of SSEGeometries corresponding to SSEs in this domain
      virtual util::ShPtrVector< SSEGeometry> GetSSEGeometries() const;

      //! @brief create and return SSEGeometries for all SSEs in this Domain
      //! @param SS_TYPE specific SSTYPE
      //! @return ShPtrVector of SSEGeometries corresponding to SSEs in this domain
      virtual util::ShPtrVector< SSEGeometry> GetSSEGeometries( const biol::SSType &SS_TYPE) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief checks if domain already contains THIS_SSE
      //! @param THIS_SSE SSE of interest
      //! @return whether domain already contains THIS_SSE
      virtual bool DoesContain( const SSE &THIS_SSE) const;

      //! @brief returns true if the domain has no SSEs in it
      //! @return true if the domain has no SSEs in it
      virtual bool IsEmpty() const;

      //! @brief checks if domain already contains this THIS_SSE or any overlapping SSE with THIS_SSE
      //! @param THIS_SSE SSE to be searched for
      //! @return if domain already contains this THIS_SSE or any overlapping SSE with THIS_SSE
      virtual bool DoesContainOverlapping( const SSE &THIS_SSE) const;

      //! @brief returns SiPtrList of sses that have short loops ( at most MAX_LOOP_LENGTH) to provided TARGET_SSE
      //! @param TARGET_SSE SSE for which short loop connecting SSEs are being search
      //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
      //! @return SiPtrList of sses that have short loops to provided TARGET_SSE
      virtual util::SiPtrList< const SSE> GetSSEsWithShortLoops
      (
        const SSE &TARGET_SSE,
        const size_t MAX_LOOP_LENGTH
      ) const;

      //! @brief selects from provided SSE_LIST, sses that have short loops ( <=MAX_LOOP_LENGTH) to SSEs in Domain
      //! @param SSE_LIST list of SSEs on which the selection will be done
      //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
      //! @return Subset of SSE_LIST that has short loops to SSEs in this domain
      virtual util::SiPtrList< const SSE> GetSSEsWithShortLoops
      (
        const util::SiPtrList< const SSE> &SSE_LIST,
        const size_t MAX_LOOP_LENGTH
      ) const;

      //! @brief returns pairs of sses that have short loops (at most MAX_LOOP_LENGTH) between each other
      //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
      //! @return list of pairs of sses that have short loops (at most MAX_LOOP_LENGTH) between
      //! @return each other
      storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > > GetSSEsWithShortLoops
      (
        const size_t MAX_LOOP_LENGTH
      ) const;

      //! @brief returns the SSE before and SSE after given TARGET_SSE in this Domain
      //! @param TARGET_SSE SSE of interest
      //! @return the SSE before and SSE after given TARGET_SSE in this Domain
      storage::VectorND< 2, util::SiPtr< const SSE> > GetNeighborSSEs
      (
        const SSE &TARGET_SSE
      ) const;

      //! @brief get the directly adjacent sses by seqid
      //! @param TARGET_SSE SSE of interest
      //! @return the SSE before and SSE after given TARGET_SSE in this Domain by seqid
      storage::VectorND< 2, util::SiPtr< const SSE> > GetAdjacentSSEs
      (
        const SSE &TARGET_SSE
      ) const;

      //! @brief returns SiPtrList of SSEs that overlap with the given TARGET_SSE
      //! @param TARGET_SSE SSE for which overlapping SSEs are searched for
      //! @return SiPtrList of SSEs that overlap with the given TARGET_SSE
      util::SiPtrList< const SSE> GetOverlappingSSEs
      (
        const SSE &TARGET_SSE
      ) const;

      //! @brief get Identification of this Domain
      //! @return string with identification
      std::string GetIdentification() const;

    }; // class DomainInterface

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_DOMAIN_INTERFACE_H_
