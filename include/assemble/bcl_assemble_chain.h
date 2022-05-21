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

#ifndef BCL_ASSEMBLE_CHAIN_H_
#define BCL_ASSEMBLE_CHAIN_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_domain_interface.h"
#include "bcl_assemble_sse.h"
#include "bcl_assemble_sse_compare.h"
#include "coord/bcl_coord_movable_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Chain
    //! @brief class Chain for representing biological unit chain, stores list of SSEs and the full sequence
    //! @details Chain class is derived from Domain ( set of SSEs) and also has the full sequence. It main functionality
    //! is keeping the full sequence and the SSEs associated since it does not encapsulate any other significant
    //! functionality in addition to ones from Domain
    //!
    //! @see @link example_assemble_chain.cpp @endlink
    //! @author woetzen, karakam
    //! @date Nov 5, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Chain :
      public DomainInterface
    {

    private:

    //////////
    // data //
    //////////

      //! sequence of the chain
      util::ShPtr< biol::AASequence> m_Sequence;

      //! Set of ShPtr to SSEs that form domain
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> m_Data;

    public:

      //! @typedef iterator over sses in chain
      typedef storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator const_ierator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Chain();

      //! @brief construct from util::ShPtr to SEQUENCE and ShPtrVector of SSEs
      //! @param SP_SEQUENCE util::ShPtr to sequence
      //! @param SSE_VECTOR ShPtrVector of SSEs
      Chain( const util::ShPtr< biol::AASequence> &SP_SEQUENCE, const util::ShPtrVector< SSE> &SSE_VECTOR);

      //! @brief  construct from util::ShPtr to SEQUENCE and DOMAIN
      //! @param SP_SEQUENCE util::ShPtr to sequence
      //! @param SOURCE_DOMAIN Domain which contains SSEs
      Chain( const util::ShPtr< biol::AASequence> &SP_SEQUENCE, const Domain &SOURCE_DOMAIN);

      //! @brief construct from util::ShPtr to SEQUENCE with no SSE information
      //! @param SP_SEQUENCE util::ShPtr to sequence
      Chain( const util::ShPtr< biol::AASequence> &SP_SEQUENCE);

      //! @brief copy constructor
      //! @param CHAIN_RHS Chain to be copied
      Chain( const Chain &CHAIN_RHS);

      //! @brief virtual copy constructor
      //! @return pointer to a new copy of this Chain
      Chain *Clone() const;

      //! @brief hardcopy
      Chain *HardCopy() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @brief the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the identification of the sequence and sses
      //! @return the identification of the sequence and sses
      std::string GetIdentification() const;

      //! @brief returns ChainID
      //! @return ChainID
      char const &GetChainID() const
      {
        return m_Sequence->GetChainID();
      }

      //! @brief set chainID to CHAINID
      //! @param CHAINID chainID to be set to
      void SetChainID( const char CHAINID);

      //! @brief return const reference to ShPtr to sequence
      //! @return const reference to ShPtr to sequence
      util::ShPtr< biol::AASequence> const &GetSequence() const
      {
        return m_Sequence;
      }

      //! @brief return non-const reference to ShPtr to sequence
      //! @return non-const reference to ShPtr to sequence
      util::ShPtr< biol::AASequence> &GetSequence()
      {
        return m_Sequence;
      }

      //! @brief return const m_Data
      //! @return const m_Data
      const storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> &GetData() const
      {
        return m_Data;
      }

      //! @brief get SiPtrVector of SSEs
      //! @return SiPtrVector of SSEs
      util::SiPtrVector< const SSE> GetSSEs() const;

      //! @brief returns all SSEs in domain of given SSTYPE in a util::SiPtrVector
      //! @param SS_TYPE specific SSTYPE
      //! @return all SSEs in domain of given SSTYPE in a util::SiPtrVector
      util::SiPtrVector< const SSE> GetSSEs( const biol::SSType &SS_TYPE) const;

      //! @brief returns SiPtrVector of all SSEs of given SS_TYPES in the set
      //! @param SS_TYPES set of SSTypes of interest
      //! @brief returns SiPtrVector of all SSEs of given SS_TYPES in the set
      util::SiPtrVector< const SSE> GetSSEs( const storage::Set< biol::SSType> &SS_TYPES) const;

      //! @brief find and to return the ShPtr for the given SSE
      //! @param SSE_TO_SEARCH SSE of interest
      //! @return ShPtr to corresponding SSE, otherwise an empty ShPtr
      const util::ShPtr< SSE> &FindSSE( const SSE &SSE_TO_SEARCH) const;

      //! @brief return total number of sses
      //! @return total number of sses
      size_t GetNumberSSEs() const
      {
        return m_Data.GetSize();
      }

      //! @brief return number of SSE of specified SSTYPE
      //! @param SS_TYPE specific SSTYPE
      //! @return number of SSE of specified SSTYPE
      size_t GetNumberSSE( const biol::SSType &SS_TYPE) const;

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      linal::Vector3D GetCenter() const;

      //! @brief returns the number of amino acids in the chain
      //! @return the number of amino acids in the chain
      size_t GetNumberAAs() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief translates the coordinates of all SSEs by the supplied translation vector
      //! @param TRANSLATION_VECTOR_3D Translation vector to be applied
      void Translate( const linal::Vector3D &TRANSLATION_VECTOR_3D);

      //! @brief transforms the coordinates of all SSEs according to given transformation matrix
      //! @param TRANSFORMATION_MATRIX_3D transformation matrix to be applied
      void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D);

      //! @brief rotate the SSE by a given rotation matrix
      //! @param ROTATION_MATRIX_3D rotation matrix to be applied
      void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D);

      //! @brief insert the given NEW_SSE into this chain
      //! @param NEW_SSE SSE to be inserted
      //! @return whether insertion was successful
      bool Insert( const util::ShPtr< SSE> &NEW_SSE);

      //! @brief insert the SSEs in the given NEW_SSE_VECTOR into this chain
      //! @param NEW_SSE_VECTOR Vector of SSEs to be inserted
      //! @return whether insertion was successful
      bool Insert( const util::ShPtrVector< SSE> &NEW_SSE_VECTOR);

      //! @brief insert the SSEs from the given domain
      //! @param NEW_DOMAIN Domain from which SSEs should be inserted
      //! @return whether insertion was successful
      bool Insert( const Domain &NEW_DOMAIN);

      //! @brief replace the given SP_SSE with already existing one
      //! @param SP_SSE ShPtr pointing to the SSE to be replaced
      bool Replace( const util::ShPtr< SSE> &SP_SSE);

      //! @brief replace all SSEs that overlap with SP_SSE with SP_SSE
      //! @param SP_SSE ShPtr to SSE to be inserted
      //! @return whether replacement succeeded
      bool ReplaceWithOverlapping( const util::ShPtr< SSE> &SP_SSE);

      //! @brief remove given SSELEMENT from the domain
      //! @param SSELEMENT SSE to be removed
      bool Remove( const SSE &SSELEMENT);

      //! @brief sets positions of all SSEs to ideal conformation w/wo superimposing with prior coordinates
      //! @param KEEP_POSITION flag to indicate whether to original body information of SSEs should be reserved
      void SetToIdealConformation( const bool KEEP_POSITION = true);

      //! @brief chop all SSE elements of that model in pieces of the sizes defined by MIN_SSE_LENGTHS
      //! @param MIN_SSE_LENGTHS VectorND of sizes that defined min size for each SSType
      void ChopSSEs( const storage::VectorND< 3, size_t> &MIN_SSE_LENGTHS);

      //! @brief checks if domain already contains THIS_SSE
      //! @param THIS_SSE SSE of interest
      //! @return whether domain already contains THIS_SSE
      bool DoesContain( const SSE &THIS_SSE) const;

      //! @brief returns true if the domain has no SSEs in it
      //! @return true if the domain has no SSEs in it
      bool IsEmpty() const
      {
        return m_Data.IsEmpty();
      }

      //! @brief checks if domain already contains this THIS_SSE or any overlapping SSE with THIS_SSE
      //! @param THIS_SSE SSE to be searched for
      //! @return if domain already contains this THIS_SSE or any overlapping SSE with THIS_SSE
      bool DoesContainOverlapping( const SSE &THIS_SSE) const;

      //! @brief AddLoops generates loop SSEs for the Chain and assigns them undefined coordinates
      //! @param UNDEFINED_COORDINATES create loop with undefined coordinates, or with coordinates form the member sequence
      //! @param MERGE_CONSECUTIVE_SSES merge consecutive sses of given type
      //! @param SS_TYPE the sstype of the sses to be merged
      void AddLoops
      (
        const bool UNDEFINED_COORDINATES,
        const bool MERGE_CONSECUTIVE_SSES,
        const biol::SSType &SS_TYPE = biol::GetSSTypes().COIL
      );

      //! @brief join following ( progressing sequence id) SSEs of given SS_TYPE into one SSE
      //! @param SS_TYPE SSType of interest
      //! @param TEST_PEPTIDE_BOND join only, if SSEs are connected by peptide bond
      void Join( const biol::SSType &SS_TYPE, const bool TEST_PEPTIDE_BOND);

      //! @brief filters the current chain by given minimum SSE sizes
      //! @param MIN_SSE_SIZES minimum SSE sizes to filter the chain by
      void FilterByMinSSESizes( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES);

      //! @brief let the aadata of the each SSE point to the corresponding data in the chain - determined by pdbID
      void ConnectSSEToChainData();

      //! @brief Replace SSEs with those drawn from the pool
      void AdoptSSEsMaintainCoordinates( const util::SiPtrVector< const SSE> &SSES);

      //! @brief Get SSE hash string to aid in identifying similar chains
      std::string GetSSEHashString() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief equal operator
      //! @param CHAIN_RHS Chain to be assigned
      //! @return this chain after being assigned to given CHAIN_RHS
      Chain &operator =( const Chain &CHAIN_RHS);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read Chain from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write Chain to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class Chain

    //! @brief construct a chain from a sequence with SSEs
    //! the phi and psi angles of the backbone are calculated and depending on the position in the ramachandran plot,
    //! SSEs are constructed
    //! @author woetzen
    //! @param SP_SEQUENCE ShPtr to sequence- at least the backbone atoms, otherwise they will be identified as loops
    //! @return Chain with given sequence and identified sses
    BCL_API Chain ConstructChainWithSSEsFromConformation( const util::ShPtr< biol::AASequence> &SP_SEQUENCE);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ChainLessThan
    //! @brief has operator for returning whether one chain is less than another
    //!
    //! @remarks example unnecessary
    //! @author weinerbe
    //! @date Dec 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct BCL_API ChainLessThan
    {
      //! @brief return whether one chain  is less than another
      //! @param CHAIN_LHS first chain
      //! @param CHAIN_RHS second chain
      //! @return whether one chain is less than another
      bool operator()( const Chain &CHAIN_LHS, const Chain &CHAIN_RHS) const;

      //! @brief return whether one chain is less than another
      //! @param PTR_CHAIN_LHS first chain
      //! @param PTR_CHAIN_RHS second chain
      //! @return whether one chain is less than another
      bool operator()
      (
        const util::PtrInterface< Chain> &PTR_CHAIN_LHS,
        const util::PtrInterface< Chain> &PTR_CHAIN_RHS
      ) const;

    }; // struct ChainLessThan

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_CHAIN_H_
