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

#ifndef BCL_ASSEMBLE_PROTEIN_MODEL_H_
#define BCL_ASSEMBLE_PROTEIN_MODEL_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_chain.h"
#include "bcl_assemble_protein_model_data.h"
#include "signal/bcl_signal.h"
#include "util/bcl_util_si_ptr_list.h"

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModel
    //! @brief class that represents a biological protein structure as a model within the bcl.
    //! @details The idea of a model bases on the reduced geometrical representation of a protein by its ideal secondary structure
    //! elements (sse, cylinder for helix, box for strand) that are assembled in three dimensional space.
    //! Depending on the resolution (complexity) of the energy function that determines the energy of the model, those
    //! sses will have aminoacids of different complexity (atoms, backbone atoms, side chains).
    //! The protein model just provides the collected functionality like moving the entity-collection of sses, replacing
    //! sses, getting all atoms in the model, and it organizes the sses - depending on their sequence - in different
    //! chains according to their chain id.
    //! For the Monte Carlo procedure it also is smart enough to make only copies of pointers of chains, if the chain is
    //! changed. This will ensure that identical data is not duplicated, and scoring only needs to rescore changed
    //! interactions.
    //!
    //! @see @link example_assemble_protein_model.cpp @endlink
    //! @author woetzen, karakam, meilerj, staritrd, fischea
    //! @date 23.04.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModel :
      public DomainInterface
    {

    protected:

    //////////
    // data //
    //////////

      //! this is the complete Sequence belonging to this model
      util::ShPtrVector< Chain> m_Chains;

      //! this is the data associated with the protein model
      util::ShPtr< ProteinModelData> m_ProteinModelData;

      //! central amino acids of fragments
      storage::Vector< storage::Pair< int, int> > m_FragmentIDs;

      //! signal handler for destructor for SSE
      mutable signal::Signal1< const ProteinModel &> m_ChangeSignal;

      //! signal handler for destructor for SSE
      mutable signal::Signal1< const ProteinModel &> m_DestructorSignal;

    public:

    ///////////
    // types //
    ///////////

      //! @typedef iterator for chains
      typedef util::ShPtrVector< Chain>::const_iterator const_iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinModel();

      //! @brief construct from ShPtrVector of Chains CHAINS
      //! @param CHAINS ShPtrVector of Chains
      ProteinModel( const util::ShPtrVector< Chain> &CHAINS);

      //! @brief construct from util::ShPtr< Chain> CHAIN
      //! @param CHAIN ShPtr to Chain
      ProteinModel( const util::ShPtr< Chain> &CHAIN);

      //! @brief copy constructor
      //! @return pointer a new ProteinModel copied from this model
      virtual ProteinModel *Clone() const;

      //! @brief hard copy constructor
      //! @return a ProteinModel with chains hard copied from that model
      virtual ProteinModel *HardCopy() const;

      //! @brief empty copy constructor
      //! @return a ProteinModel that is empty
      virtual ProteinModel *Empty() const;

      //! @brief destructor
      virtual ~ProteinModel();

      //! @brief hard copy a protein model by hard copying its chains
      static ProteinModel HardCopy( const ProteinModel &PROTEIN_MODEL);

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the identification of the chains
      //! @return the identification of the chains
      std::string GetIdentification() const;

      //! @brief access the protein model data object which contains miscellaneous data for that protein model
      //! @return SHPtr to the protein model data object
      virtual const util::ShPtr< ProteinModelData> &GetProteinModelData() const;

      //! @brief set a new protein model data object
      //! @param SP_PROTEIN_MODEL_DATA ShPtr to a new protein model data object
      virtual void SetProteinModelData( const util::ShPtr< ProteinModelData> &SP_PROTEIN_MODEL_DATA);

      //! @brief sets the sse pool within the protein model data
      //! @param SSE_POOL the sse pool that will be set in the protein model data
      virtual void SetSSEPoolData( const util::ShPtr< SSEPool> &SSE_POOL);

      //! @brief returns number of structure elements in model
      //! @return number of SSEs in model
      virtual size_t GetNumberOfChains() const
      {
        return m_Chains.GetSize();
      }

      //! @brief returns total Number of SSEs
      //! @return total Number of SSEs
      virtual size_t GetNumberSSEs() const;

      //! @brief return number of SSE of specified SSTYPE
      //! @param SS_TYPE specific SSTYPE
      //! @return number of SSE of specified SSTYPE
      virtual size_t GetNumberSSE( const biol::SSType &SS_TYPE) const;

      //! @brief returns const reference to ShPtrVector of all chains in this protein model
      //! @return const reference to ShPtrVector of all chains in this protein model
      virtual const util::ShPtrVector< Chain> &GetChains() const
      {
        return m_Chains;
      }

      //! @brief returns non-const reference to ShPtrVector of all chains in this protein model
      //! @return non-const reference to ShPtrVector of all chains in this protein model
      virtual util::ShPtrVector< Chain> &GetChains()
      {
        return m_Chains;
      }

      //! @brief get specific chains
      //! @param CHAIN_IDS specific chain IDs as a string
      //! @return specific chains requested
      virtual util::ShPtrVector< Chain> GetChains( const std::string &CHAIN_IDS) const;

      //! @brief returns chains without SSE information
      //! @return chains without SSE information
      virtual util::ShPtrVector< Chain> GetEmptyChains() const;

      //! @brief get all chain ids
      //! @return vector with all chain ids in protein model
      virtual std::string GetChainIDs() const;

      //! @brief returns const ShPtr to Chain with the supplied chain id
      //! @param CHAIN_ID Chain ID of the chain that is being searched
      //! @return const ShPtr to Chain with the supplied chain id
      virtual const util::ShPtr< Chain> &GetChain( const char CHAIN_ID) const;

      //! @brief returns const ShPtr to Chain with the supplied chain id
      //! @param CHAIN_ID Chain ID of the chain that is being searched
      //! @return const ShPtr to Chain with the supplied chain id
      virtual util::ShPtr< Chain> &GetChain( const char CHAIN_ID);

      //! @brief returns SiPtrVector of AASequences of all chains
      //! @return SiPtrVector of AASequences of all chains
      virtual util::SiPtrVector< const biol::AASequence> GetSequences() const;

      //! @brief returns SiPtrVector of all SSEs of given SS_TYPE
      //! @param SS_TYPE SSType of interest
      //! @return SiPtrVector of all SSEs of given SS_TYPE
      virtual util::SiPtrVector< const SSE> GetSSEs( const biol::SSType &SS_TYPE) const;

      //! @brief returns SiPtrVector of all SSEs of given SS_TYPES in the set
      //! @param SS_TYPES set of SSTypes of interest
      //! @brief returns SiPtrVector of all SSEs of given SS_TYPES in the set
      virtual util::SiPtrVector< const SSE> GetSSEs( const storage::Set< biol::SSType> &SS_TYPES) const;

      //! @brief returns SiPtrVector of all SSEs in all chains
      //! @return SiPtrVector of all SSEs in all chains
      virtual util::SiPtrVector< const SSE> GetSSEs() const;

      //! @brief returns the SSE containing the given residue
      //! @param RESIDUE residue for which to return the SSE
      //! @return SSE containing the given residue
      virtual util::SiPtr< const SSE> GetSSE( const biol::AABase &RESIDUE) const;

      //! @brief returns ShPtrVector of all SSE geometries in all chainssp_protein_model
      //! @return ShPtrVector of all SSE geometries in all chains
      virtual util::ShPtrVector< SSEGeometry> GetSSEGeometries() const;

      //! @brief returns ShPtrVector of all SSE fragments
      //! @return ShPtrVector of all SSE fragments
      virtual util::ShPtrVector< SSEGeometryInterface> GetFragments() const;

      //! @brief returns all SSEs in Chains in a single Domain
      //! @return all SSEs in Chains in a single Domain
      virtual Domain GetSSEsAsDomain() const;

      //! @brief returns all SSEs of SS_TYPE in Chains in a single Domain
      //! @param SS_TYPE SSType of interest
      //! @return all SSEs of specified SS_TYPE in Chains in a single Domain
      virtual Domain GetSSEsAsDomain( const biol::SSType &SS_TYPE) const;

      //! @brief find and to return the ShPtr for the given SSE
      //! @param SSE_TO_SEARCH SSE of interest
      //! @return ShPtr to corresponding SSE, otherwise an empty ShPtr
      virtual const util::ShPtr< SSE> &FindSSE( const SSE &SSE_TO_SEARCH) const;

      //! @brief returns concatenated vector of all amino acids in all SSEs in all chains(excludes amino acids in loops)
      //! @return concatenated vector of all amino acids in all SSEs in all chains(excludes amino acids in loops)
      virtual util::SiPtrVector< const biol::AABase> GetAminoAcids() const;

      //! @brief returns all atoms in a SiPtrVector
      //! @return all atoms in a SiPtrVector
      virtual util::SiPtrVector< const biol::Atom> GetAtoms() const;

      //! @brief returns all atoms of specified types in a SiPtrVector
      //! @param ATOM_TYPES Set of AtomTypes of interest
      //! @return all atoms of specified types in a SiPtrVector
      virtual util::SiPtrVector< const biol::Atom> GetAtoms
      (
        const storage::Set< biol::AtomType> &ATOM_TYPES
      ) const;

      //! @brief returns coordinates for all atoms in a SiPtrVector
      //! @return coordinates for all atoms in a SiPtrVector
      virtual util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates() const;

      //! @brief returns coordinates for all atoms of specified types in a SiPtrVector
      //! @param ATOM_TYPES Set of AtomTypes of interest
      //! @return coordinates for all atoms of specified types in a SiPtrVector
      virtual util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates
      (
        const storage::Set< biol::AtomType> &ATOM_TYPES
      ) const;

      //! @brief returns defined coordinates for all atoms of specified types in a SiPtrVector
      //! @param ATOM_TYPES Set of AtomTypes of interest
      //! @return defined coordinates for all atoms of specified types in a SiPtrVector
      virtual util::SiPtrVector< const linal::Vector3D> GetDefinedAtomCoordinates
      (
        const storage::Set< biol::AtomType> &ATOM_TYPES
      ) const;

      //! @brief returns the center of the protein model
      //! @return the center of the protein model
      virtual linal::Vector3D GetCenter() const;

      //! @brief return centers of SSES
      //! @return center of mass of centers of sses
      virtual linal::Vector3D GetCenterOfSSEs() const;

      //! @brief return the center of mass of this model
      //! @return the center of mass of this model
      virtual linal::Vector3D GetCenterOfMass() const
      {
        // center of mass of all atoms coordinates, skipping undefined coordinates
        return coord::CenterOfMass( GetAtomCoordinates(), true);
      }

      //! @brief access to the change signal handler
      //! @return const ref to the SignalHandler that emits on changes
      signal::Signal1< const ProteinModel &> &GetChangeSignal() const
      {
        return m_ChangeSignal;
      }

      //! @brief access to the destructor signal handler
      //! @return const ref to the SignalHandler that emits on destruction
      signal::Signal1< const ProteinModel &> &GetDestructorSignal() const
      {
        return m_DestructorSignal;
      }

      //! @brief gives the conformations that are in this protein
      //! @return protein ensemble which is the conformations making up this protein
      virtual const ProteinEnsemble &GetConformationalEnsemble() const;

      //! @brief sets the conformations that are in this protein
      //! @param ENSEMBLE protein ensemble which is the conformations making up this protein
      virtual void SetConformationalEnsemble( const ProteinEnsemble &ENSEMBLE);

      //! @brief gives the ids of SSE fragments
      //! @return vector of IDs of SSE fragments
      storage::Vector< storage::Pair< int, int> > GetFragmentIDs() const
      {
        return m_FragmentIDs;
      }

      //! @brief gives the number of amino acids in the protein
      //! @return the number of amino acids in the protein
      size_t GetNumberAAs() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief transforms the coordinates of the entire model with given TransformationMatrix3D
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      virtual void Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D);

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION vector to be used in translation
      virtual void Translate( const linal::Vector3D &TRANSLATION);

      //! @brief rotate the object by a given RotationMatrix3D
      //! @param ROTATION_MATRIX_3D RotationMatrix3D to be used in rotation
      virtual void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D);

      //! @brief replaces SSEs in the protein model with the given SSE
      //! @detail SSEs in the protein model that are overlapping with the given SSE are either completely
      //! replaced or resized in order to not overlap with the given SSE
      //! @param SP_SSE SSE that is to be inserted into the protein model
      virtual void ReplaceResize( const util::ShPtr< SSE> &SP_SSE);

      //! @brief adds a chain to the model
      //! @param SP_CHAIN ShPtr to Chain to be inserted
      virtual void Insert( const util::ShPtr< Chain> &SP_CHAIN);

      //! @brief adds a structure element to the model
      //! @brief SP_SSE ShPtr to the SSE to be inserted
      //! @return whether insertion succeeded
      virtual bool Insert( const util::ShPtr< SSE> &SP_SSE);

      //! @brief replaces given SSE if it was in one of the chains
      //! @param SP_SSE ShPtr to the SSE of interest
      //! @return whether replacement succeeded
      virtual bool Replace( const util::ShPtr< SSE> &SP_SSE);

      //! @brief replaces given chain if it was in the protein model
      //! @param SP_CHAIN ShPtr to the chain of interest
      //! @return whether replacement succeeded
      virtual bool Replace( const util::ShPtr< Chain> &SP_CHAIN);

      //! @brief replaces vector of SSEs in the chains
      //! @brief SSELEMENTS ShPtrVector of SSEs to be replaced
      //! @return whether replacement succeeded
      virtual bool Replace( const util::ShPtrVector< SSE> &SSELEMENTS);

      //! @brief replace all SSEs that overlap with SP_SSE with SP_SSE
      //! @param SP_SSE ShPtr to SSE to be inserted
      //! @return whether replacement succeeded
      virtual bool ReplaceWithOverlapping( const util::ShPtr< SSE> &SP_SSE);

      //! @brief replaces given set of SSEs in the chains
      //! @param SSELEMENTS set of SSEs to be replaced
      //! @return whether replacement succeeded
      virtual bool Replace( const storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> &SSELEMENTS);

      //! @brief remove given SSE
      //! @brief THIS_SSE SSE of interest
      //! @return whether removal succeeded
      virtual bool Remove( const SSE &THIS_SSE);

      //! @brief sets all positions of all SSEs to ideal conformation and superimpose them with prior coordinates
      //! depending on flag KEEP_POSITION
      //! @param KEEP_POSITION flag to determine superimposition with prior coordinates
      virtual void SetToIdealConformation( const bool KEEPPOSITION = true);

      //! @brief chop all sselements of that model in pieces of the given size vector MIN_SSE_LENGTHS
      //! @param MIN_SSE_LENGTHS minimal sizes for chopping for each SSType
      virtual void ChopSSEs( const storage::VectorND< 3, size_t> &MIN_SSE_LENGTHS);

      //! @brief checks whether the given SSE is already in the model or is overlapping with one
      //! @param THIS_SSE SSE of interest
      //! @return whether the given SSE is already in the model or is overlapping with one
      virtual bool DoesContain( const SSE &THIS_SSE) const;

      //! @brief returns pairs of sses that have short loops (at most MAX_LOOP_LENGTH) between each other
      //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
      //! @return list of pairs of sses that have short loops (at most MAX_LOOP_LENGTH) between
      //! @return each other
      virtual storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > > GetSSEsWithShortLoops
      (
        const size_t MAX_LOOP_LENGTH
      ) const;

      //! @brief returns SiPtrList of sses that have short loops ( at most MAX_LOOP_LENGTH) to provided TARGET_SSE
      //! @param TARGET_SSE SSE for which short loop connecting SSEs are being search
      //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
      //! @return SiPtrList of sses that have short loops to provided TARGET_SSE
      virtual util::SiPtrList< const SSE> GetSSEsWithShortLoops
      (
        const SSE &TARGET_SSE,
        const size_t MAX_LOOP_LENGTH
      ) const;

      //! @brief selects from provided SSE_LIST, sses that have short loops ( <=MAX_LOOP_LENGTH) to SSEs in ProteinModel
      //! @param SSE_LIST list of SSEs on which the selection will be done
      //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
      //! @return Subset of SSE_LIST that has short loops to SSEs in this ProteinModel
      virtual util::SiPtrList< const SSE> GetSSEsWithShortLoops
      (
        const util::SiPtrList< const SSE> &SSE_LIST,
        const size_t MAX_LOOP_LENGTH
      ) const;

      //! @brief returns the SSE before and SSE after given TARGET_SSE in this ProteinModel
      //! @param TARGET_SSE SSE of interest
      //! @return the SSE before and SSE after given TARGET_SSE in this ProteinModel
      virtual storage::VectorND< 2, util::SiPtr< const SSE> > GetNeighborSSEs
      (
        const SSE &TARGET_SSE
      ) const;

      //! @brief returns SSEs in this ProteinModel which overlap with given TARGET_SSE
      //! @brief TARGET_SSE SSE for which overlaps are going to be searched
      //! @return SSEs in this ProteinModel which overlap with given TARGET_SSE
      virtual util::SiPtrList< const SSE> GetOverlappingSSEs
      (
        const SSE &TARGET_SSE
      ) const;

      //! @brief returns SSEs in this ProteinModel which overlap with one or more SSEs in the given SSE_LIST
      //! @brief SSE_LIST SiPtrList of SSEs for which overlaps are going to be searched
      //! @return SSEs in this ProteinModel which overlap with one or more SSEs in the given SSE_LIST
      virtual util::SiPtrList< const SSE> GetOverlappingSSEs
      (
        const util::SiPtrList< const SSE> &SSE_LIST
      ) const;

      //! @brief AddLoops generates loop SSEs for the protein model and assigns them zero coordinates
      //! @param UNDEFINED_COORDINATES create loop with undefined coordinates, or with coordinates form the member sequence
      //! @param MERGE_CONSECUTIVE_SSES merge consecutive sses of given type
      //! @param SS_TYPE hte sstype of the sses to be merged
      virtual void AddLoops
      (
        const bool UNDEFINED_COORDINATES,
        const bool MERGE_CONSECUTIVE_SSES,
        const biol::SSType &SS_TYPE = biol::GetSSTypes().COIL
      );

      //! @brief filters the current protein model by given minimum SSE sizes
      //! @param MIN_SSE_SIZES minimum SSE sizes to filter the model by
      virtual void FilterByMinSSESizes
      (
        const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES
      );

      //! @brief join following ( progressing sequence id) SSEs of given SS_TYPE into one SSE
      //! @param SS_TYPE SSType of interest
      //! @param TEST_PEPTIDE_BOND join only, if SSEs are connected by peptide bond
      virtual void Join( const biol::SSType &SS_TYPE, const bool TEST_PEPTIDE_BOND);

      //! @brief let the aadata of the each sse point to the corresponding data in the chain - determined by pdbID
      virtual void ConnectSSEToChainData();

      //! @brief Replace SSEs with those drawn from the pool
      virtual void AdoptSSEsMaintainCoordinates( const util::SiPtrVector< const SSE> &SSES);

      //! @brief Get SSE hash string to aid in identifying similar proteins
      std::string GetSSEHashString() const;

      //! @brief calculates all pairwise interactions of amino acids within SSE fragments
      //! @param SIZE number of amino acids in the chain
      //! @return matrix of pairwise interactions of amino acids within SSE fragments
      const linal::Matrix< double> CalculateFragmentInteractions() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator
      //! @param PROTEIN_MODEL ProteinModel to be copied
      //! @return This model after all members are assigned to values from PROTEIN_MODEL
      virtual ProteinModel &operator =( const ProteinModel &PROTEIN_MODEL);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read ProteinModel from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write ProteinModel to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ProteinModel

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_PROTEIN_MODEL_H_
