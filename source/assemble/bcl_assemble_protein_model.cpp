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
#include "assemble/bcl_assemble_protein_model.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_topology_combined.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_sse_geometry_packing_pickers.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "storage/bcl_storage_row.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    //! @class TransformationsCompare
    //! @brief compares transformations
    struct BCL_API TransformationsCompare
    {
      //! @brief compares transformations
      //! @param TRIPLET_A first transformation
      //! @param TRIPLET_B second transformation
      //! @return comparison
      bool operator()
      (
        const storage::Triplet< char, char, math::TransformationMatrix3D> &TRIPLET_A,
        const storage::Triplet< char, char, math::TransformationMatrix3D> &TRIPLET_B
      )
      {
        return TRIPLET_A.Second() < TRIPLET_B.Second();
      }
    };

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ProteinModel::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModel())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModel::ProteinModel() :
      m_Chains(),
      m_ProteinModelData( new ProteinModelData())
    {
    }

    //! @brief construct from ShPtrVector of Chains CHAINS
    //! @param CHAINS ShPtrVector of Chains
    ProteinModel::ProteinModel( const util::ShPtrVector< Chain> &CHAINS) :
      m_Chains( CHAINS),
      m_ProteinModelData( new ProteinModelData())
    {
    }

    //! @brief construct from util::ShPtr< Chain> CHAIN
    //! @param CHAIN ShPtr to Chain
    ProteinModel::ProteinModel( const util::ShPtr< Chain> &CHAIN) :
      m_Chains( 1, CHAIN),
      m_ProteinModelData( new ProteinModelData())
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer a new ProteinModel copied from this model
    ProteinModel *ProteinModel::Clone() const
    {
      return new ProteinModel( *this);
    }

    //! @brief hard copy constructor
    //! @return a ProteinModel with chains hard copied from that model
    ProteinModel *ProteinModel::HardCopy() const
    {
      return new ProteinModel( HardCopy( *this));
    }

    //! @brief hard copy a protein model by hard copying its chains
    ProteinModel ProteinModel::HardCopy( const ProteinModel &PROTEIN_MODEL)
    {
      // new protein model
      ProteinModel hard_copy;
      hard_copy.SetProteinModelData( PROTEIN_MODEL.m_ProteinModelData->HardCopy());

      // iterate over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.m_Chains.Begin()), chain_itr_end( PROTEIN_MODEL.m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // push back a shptr to a hardcopied chain
        hard_copy.m_Chains.PushBack( util::ShPtr< Chain>( ( *chain_itr)->HardCopy()));
      }

      // end
      return hard_copy;
    }

    //! @brief empty copy constructor
    //! @return a ProteinModel that is empty
    ProteinModel *ProteinModel::Empty() const
    {
      return new ProteinModel();
    }

    //! @brief destructor
    ProteinModel::~ProteinModel()
    {
      m_DestructorSignal.Emit( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModel::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the identification of the chains
    //! @return the identification of the chains
    std::string ProteinModel::GetIdentification() const
    {
      // initialize string
      std::string identification;

      // iterate through the chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // add the chain identification to the string
        identification += ( *chain_itr)->GetIdentification();
      }

      // end
      return identification;
    }

    //! @brief access the protein model data object which contains miscellaneous data for that protein model
    //! @return SHPtr to the protein model data object
    const util::ShPtr< ProteinModelData> &ProteinModel::GetProteinModelData() const
    {
      return m_ProteinModelData;
    }

    //! @brief set a new protein model data object
    //! @param SP_PROTEIN_MODEL_DATA ShPtr to a new protein model data object
    void ProteinModel::SetProteinModelData( const util::ShPtr< ProteinModelData> &SP_PROTEIN_MODEL_DATA)
    {
      m_ProteinModelData = SP_PROTEIN_MODEL_DATA;
      m_ChangeSignal.Emit( *this);
    }

    //! @brief sets the sse pool within the protein model data
    //! @param SSE_POOL the sse pool that will be set in the protein model data
    void ProteinModel::SetSSEPoolData( const util::ShPtr< SSEPool> &SSE_POOL)
    {
      // copy the ProteinModelData
      util::ShPtr< ProteinModelData> sp_model_data( GetProteinModelData()->HardCopy());
      if( !sp_model_data->Insert( ProteinModelData::e_Pool, SSE_POOL))
      {
        BCL_Assert
        (
          sp_model_data->Replace( ProteinModelData::e_Pool, SSE_POOL),
          "could not replace protein model data " + util::Format()( assemble::ProteinModelData::e_Pool)
        );
      }

      SetProteinModelData( sp_model_data);

      m_ChangeSignal.Emit( *this);
    }

    //! @brief returns total Number of SSEs
    //! @return total Number of SSEs
    size_t ProteinModel::GetNumberSSEs() const
    {
      // initialize counter
      size_t number_sses( 0);

      // loop over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // sum up the number of such SSEs from this chain
        number_sses += ( *chain_itr)->GetNumberSSEs();
      }
      // end
      return number_sses;
    }

    //! @brief return number of SSE of specified SSTYPE
    //! @param SS_TYPE specific SSTYPE
    //! @return number of SSE of specified SSTYPE
    size_t ProteinModel::GetNumberSSE( const biol::SSType &SS_TYPE) const
    {
      // initialize counter
      size_t number_sses( 0);

      // loop over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // sum up the number of such SSEs from this chain
        number_sses += ( *chain_itr)->GetNumberSSE( SS_TYPE);
      }
      // end
      return number_sses;
    }

    //! @brief get specific chains
    //! @param CHAIN_IDS specific chain IDs as a string
    //! @return specific chains requested
    util::ShPtrVector< Chain> ProteinModel::GetChains( const std::string &CHAIN_IDS) const
    {
      util::ShPtrVector< Chain> chains;
      for( std::string::const_iterator itr( CHAIN_IDS.begin()); itr != CHAIN_IDS.end(); ++itr)
      {
        chains.PushBack( GetChain( *itr));
      }
      return chains;
    }

    //! @brief returns chains without SSE information
    //! @return chains without SSE information
    util::ShPtrVector< Chain> ProteinModel::GetEmptyChains() const
    {
      // initialize vector of chains
      util::ShPtrVector< Chain> chains;

      // iterate over the chains in the protein model
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // pushback the chain lacking SSE information into the chains vector
        chains.PushBack( util::ShPtr< Chain>( new Chain( ( *chain_itr)->GetSequence())));
      }

      //end
      return chains;
    }

    //! @brief get all chain ids
    //! @return vector with all chain ids in protein model
    std::string ProteinModel::GetChainIDs() const
    {
      std::string chain_ids;

      //loop over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        chain_ids.push_back( ( *chain_itr)->GetChainID());
      }

      // end
      return chain_ids;
    }

    //! @brief returns const ShPtr to Chain with the supplied chain id
    //! @param CHAIN_ID Chain ID of the chain that is being searched
    //! @return const ShPtr to Chain with the supplied chain id
    const util::ShPtr< Chain> &ProteinModel::GetChain( const char CHAIN_ID) const
    {
      //loop over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // if this chain has the supplied CHAIN_ID, return it
        if( ( *chain_itr)->GetChainID() == CHAIN_ID)
        {
          return *chain_itr;
        }
      }

      // if no such chain was found message and return empty shared pointer
      BCL_MessageDbg
      (
        "No chain with the supplied chain_id '" + util::Format()( CHAIN_ID) +
        "' was found in this protein model with the chains\n" + util::Format()( GetChainIDs())
      );
      static const util::ShPtr< Chain> s_undefined_chain;
      return s_undefined_chain;
    }

    //! @brief returns const ShPtr to Chain with the supplied chain id
    //! @param CHAIN_ID Chain ID of the chain that is being searched
    //! @return const ShPtr to Chain with the supplied chain id
    util::ShPtr< Chain> &ProteinModel::GetChain( const char CHAIN_ID)
    {
      //loop over all chains
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // if this chain has the supplied CHAIN_ID, return it
        if( ( *chain_itr)->GetChainID() == CHAIN_ID)
        {
          return *chain_itr;
        }
      }

      // if no such chain was found message and return empty shared pointer
      BCL_MessageDbg
      (
        "No chain with the supplied chain_id '" + util::Format()( CHAIN_ID) +
        "' was found in this protein model with the chains\n" + util::Format()( GetChainIDs())
      );
      static util::ShPtr< Chain> s_undefined_chain;
      return s_undefined_chain;
    }

    //! @brief returns SiPtrVector of AASequences of all chains
    //! @return SiPtrVector of AASequences of all chains
    util::SiPtrVector< const biol::AASequence> ProteinModel::GetSequences() const
    {
      // initialize sequences vector
      util::SiPtrVector< const biol::AASequence> sequences;

      //loop over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        sequences.PushBack( ( *chain_itr)->GetSequence());
      }

      // end
      return sequences;
    }

    //! @brief returns SiPtrVector of all SSEs of given SS_TYPE
    //! @param SS_TYPE SSType of interest
    //! @return SiPtrVector of all SSEs of given SS_TYPE
    util::SiPtrVector< const SSE> ProteinModel::GetSSEs( const biol::SSType &SS_TYPE) const
    {
      util::SiPtrVector< const SSE> sses;

      //loop over all SSElements
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        sses.Append( ( *chain_itr)->GetSSEs( SS_TYPE));
      }

      return sses;
    }

    //! @brief returns SiPtrVector of all SSEs of given SS_TYPES in the set
    //! @param SS_TYPES set of SSTypes of interest
    //! @brief returns SiPtrVector of all SSEs of given SS_TYPES in the set
    util::SiPtrVector< const SSE> ProteinModel::GetSSEs( const storage::Set< biol::SSType> &SS_TYPES) const
    {
      util::SiPtrVector< const SSE> sses;

      //loop over all SSElements
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        sses.Append( ( *chain_itr)->GetSSEs( SS_TYPES));
      }

      return sses;
    }

    //! @brief returns SiPtrVector of all SSEs in all chains
    //! @return SiPtrVector of all SSEs in all chains
    util::SiPtrVector< const SSE> ProteinModel::GetSSEs() const
    {
      util::SiPtrVector< const SSE> sses;

      //loop over all SSElements
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        sses.Append( ( *chain_itr)->GetSSEs());
      }

      return sses;
    }

    //! @brief returns the SSE containing the given residue
    //! @param RESIDUE residue for which to return the SSE
    //! @return SSE containing the given residue
    util::SiPtr< const SSE> ProteinModel::GetSSE( const biol::AABase &RESIDUE) const
    {
      // get the chain and sequence ids
      const char chain_id( RESIDUE.GetChainID());
      const int seq_id( RESIDUE.GetSeqID());

      // find the SSE containing the given residue
      const Chain &chain( *GetChain( chain_id));
      for
      (
        auto sse_it( chain.GetData().Begin()), sse_it_end( chain.GetData().End());
        sse_it != sse_it_end;
        ++sse_it
      )
      {
        const SSE &sse( **sse_it);
        if( sse.GetFirstAA()->GetSeqID() <= seq_id && sse.GetLastAA()->GetSeqID() >= seq_id)
        {
          return *sse_it;
        }
      }

      return util::SiPtr< const SSE>();
    }

    //! @brief returns ShPtrVector of all SSE geometries in all chains
    //! @return ShPtrVector of all SSE geometries in all chains
    util::ShPtrVector< SSEGeometry> ProteinModel::GetSSEGeometries() const
    {
      // initialize vector of geometries
      util::ShPtrVector< SSEGeometry> geometry_vector;

      // iterate over the chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator itr( m_Chains.Begin()), itr_end( m_Chains.End());
        itr != itr_end; ++itr
      )
      {
        // pushback the geometries from this chain
        geometry_vector.Append( ( *itr)->GetSSEGeometries());
      }
      // end
      return geometry_vector;
    }

    //! @brief returns ShPtrVector of all SSE fragments in all chains
    //! @return ShPtrVector of all SSE fragments in all chains
    util::ShPtrVector< SSEGeometryInterface> ProteinModel::GetFragments() const
     {
       // initialize vector of geometries
       util::ShPtrVector< SSEGeometryInterface> geometry_vector;

       // iterate over SSEs and add fragments to geometry vector
       for
       (
         util::ShPtrVector< Chain>::const_iterator itr( m_Chains.Begin()), itr_end( m_Chains.End());
         itr != itr_end; ++itr
       )
       {
         util::SiPtrVector< const SSE> sses( ( *itr)->GetSSEs());
         for
         (
             util::SiPtrVector< const SSE>::const_iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
             sse_itr != sse_itr_end; sse_itr++
         )
         {
           geometry_vector.Append( ( *sse_itr)->GetFragments());
         }
       }
       // end
       return geometry_vector;
     }

    //! @brief returns all SSEs in Chains in a single Domain
    //! @return all SSEs in Chains in a single Domain
    Domain ProteinModel::GetSSEsAsDomain() const
    {
      // initialize domain to return
      Domain this_domain;

      //loop over all Chaiins
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // insert the domain
        this_domain.Insert( ( *chain_itr)->GetData());
      }

      // end
      return this_domain;
    }

    //! @brief returns all SSEs of SS_TYPE in Chains in a single Domain
    //! @param SS_TYPE SSType of interest
    //! @return all SSEs of specified SS_TYPE in Chains in a single Domain
    Domain ProteinModel::GetSSEsAsDomain( const biol::SSType &SS_TYPE) const
    {
      // initialize domain to return
      Domain this_domain;

      //loop over all Chaiins
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // iterate over all the SSEs in this chain
        for
        (
          storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          // if the type is correct
          if( ( *sse_itr)->GetType() == SS_TYPE)
          {
            // insert into domain
            this_domain.Insert( *sse_itr);
          }
        }
      }

      // end
      return this_domain;
    }

    //! @brief find and to return the ShPtr for the given SSE
    //! @param SSE_TO_SEARCH SSE of interest
    //! @return ShPtr to corresponding SSE, otherwise an empty ShPtr
    const util::ShPtr< SSE> &ProteinModel::FindSSE( const SSE &SSE_TO_SEARCH) const
    {
      // static undefined SSE ptr
      static const util::ShPtr< SSE> s_undefined_sse_ptr;

      //iterate over all chains, and find the right Chain matching sse's chain id
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        //if matching chain id
        if( ( *chain_itr)->GetChainID() == SSE_TO_SEARCH.GetChainID())
        {
          // find the SSE in the chain and return it
          return ( *chain_itr)->FindSSE( SSE_TO_SEARCH);
        }
      }

      // if this point is reached then no chain with matching chain id is found
      return s_undefined_sse_ptr;
    }

    //! @brief returns concatenated vector of all amino acids in all sses in all chains(excludes amino acids in loops)
    //! @return concatenated vector of all amino acids in all sses in all chains(excludes amino acids in loops)
    util::SiPtrVector< const biol::AABase> ProteinModel::GetAminoAcids() const
    {
      util::SiPtrVector< const biol::AABase> amino_acids;
      amino_acids.AllocateMemory( GetNumberAAs());

      //loop over all SSElements
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        amino_acids.Append( ( *chain_itr)->GetAminoAcids());
      }

      return amino_acids;
    }

    //! @brief returns all atoms in a SiPtrVector
    //! @return all atoms in a SiPtrVector
    util::SiPtrVector< const biol::Atom> ProteinModel::GetAtoms() const
    {
      util::SiPtrVector< const biol::Atom> atoms;

      //loop over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        atoms.Append( ( *chain_itr)->GetAtoms());
      }

      return atoms;
    }

    //! @brief returns all atoms of specified types in a SiPtrVector
    //! @param ATOM_TYPES Set of AtomTypes of interest
    //! @return all atoms of specified types in a SiPtrVector
    util::SiPtrVector< const biol::Atom> ProteinModel::GetAtoms
    (
      const storage::Set< biol::AtomType> &ATOM_TYPES
    ) const
    {
      // util::SiPtrVector of atoms
      util::SiPtrVector< const biol::Atom> atoms;

      // loop over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        atoms.Append( ( *chain_itr)->GetAtoms( ATOM_TYPES));
      }

      //return
      return atoms;
    }

    //! @brief returns coordinates for all atoms in a SiPtrVector
    //! @return coordinates for all atoms in a SiPtrVector
    util::SiPtrVector< const linal::Vector3D> ProteinModel::GetAtomCoordinates() const
    {
      //util::SiPtrVector of positions
      util::SiPtrVector< const linal::Vector3D> positions;

      //loop over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        positions.Append( ( *chain_itr)->GetAtomCoordinates());
      }

      //return
      return positions;
    }

    //! @brief returns coordinates for all atoms of specified types in a SiPtrVector
    //! @param ATOM_TYPES Set of AtomTypes of interest
    //! @return coordinates for all atoms of specified types in a SiPtrVector
    util::SiPtrVector< const linal::Vector3D> ProteinModel::GetAtomCoordinates
    (
      const storage::Set< biol::AtomType> &ATOM_TYPES
    ) const
    {
      //util::SiPtrVector of positions
      util::SiPtrVector< const linal::Vector3D> positions;

      //loop over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        positions.Append( ( *chain_itr)->GetAtomCoordinates( ATOM_TYPES));
      }

      //return
      return positions;
    }

    //! @brief returns defined coordinates for all atoms of specified types in a SiPtrVector
    //! @param ATOM_TYPES Set of AtomTypes of interest
    //! @return defined coordinates for all atoms of specified types in a SiPtrVector
    util::SiPtrVector< const linal::Vector3D> ProteinModel::GetDefinedAtomCoordinates
    (
      const storage::Set< biol::AtomType> &ATOM_TYPES
    ) const
    {
      // SiPtrVector that holds defined coordinates of all atoms
      util::SiPtrVector< const linal::Vector3D> defined_coordinates;

      // loop over all chains
      for
      (
        auto chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // get all atoms of ATOM_TYPE on the current chain
        util::SiPtrVector< const biol::Atom> atoms( ( *chain_itr)->GetAtoms( ATOM_TYPES));

        // loop over all atoms
        for
        (
          auto atom_itr( atoms.Begin()), atom_itr_end( atoms.End());
          atom_itr != atom_itr_end;
          ++atom_itr
        )
        {
          // only push back defined coordinates
          if( ( *atom_itr)->GetCoordinates().IsDefined())
          {
            defined_coordinates.PushBack( ( *atom_itr)->GetCoordinates());
          }
        }
      }

      // return all defined atom coordinates
      return defined_coordinates;
    }

    //! @brief returns the center of the protein model
    //! @return the center of the protein model
    linal::Vector3D ProteinModel::GetCenter() const
    {
      return coord::CenterOfMass( GetAtomCoordinates( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA)), true);
    }

    //! @brief returns the center of the centers of the sses
    //! @return center of mass of centers of sses
    linal::Vector3D ProteinModel::GetCenterOfSSEs() const
    {
      // collect center of all SSEs
      const util::SiPtrVector< const SSE> all_sses( GetSSEs());

      // initialize center
      linal::Vector3D center( 0.0);

      size_t count_sse( 0);
      // iterate over all sses
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // skip if body is undefined
        if( !( *sse_itr)->IsDefined())
        {
          continue;
        }

        center += ( *sse_itr)->GetCenter();
        ++count_sse;
      }

      //normalize
      center /= double( count_sse);

      return center;
    }

    //! @brief gives the conformations that are in this protein
    //! @return protein ensemble which is the conformations making up this protein
    const ProteinEnsemble &ProteinModel::GetConformationalEnsemble() const
    {
      static ProteinEnsemble ensemble;

      return ensemble;
    }

    //! @brief sets the conformations that are in this protein
    //! @param ENSEMBLE protein ensemble which is the conformations making up this protein
    void ProteinModel::SetConformationalEnsemble( const ProteinEnsemble &ENSEMBLE)
    {
    }

    //! @brief gives the number of amino acids in the protein
    //! @return the number of amino acids in the protein
    size_t ProteinModel::GetNumberAAs() const
    {
      // initialize size
      size_t size( 0);
      // loop over sequence to get size of protein model
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        size += ( *chain_itr)->GetNumberAAs();
      }
      return ( size);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief transforms the coordinates of the entire model with given TransformationMatrix3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void ProteinModel::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      //loop over all chains and transform them
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        ( *chain_itr)->Transform( TRANSFORMATION_MATRIX_3D);
      }

      // emit change signal
      m_ChangeSignal.Emit( *this);
    }

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION vector to be used in translation
    void ProteinModel::Translate( const linal::Vector3D &TRANSLATION)
    {
      //loop over all chains and transform them
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        ( *chain_itr)->Translate( TRANSLATION);
      }

      // emit change signal
      m_ChangeSignal.Emit( *this);
    }

    //! @brief rotate the object by a given RotationMatrix3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be used in rotation
    void ProteinModel::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      //loop over all chains and transform them
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        ( *chain_itr)->Rotate( ROTATION_MATRIX_3D);
      }

      // emit change signal
      m_ChangeSignal.Emit( *this);
    }

    //! @brief adds a chain to the model
    //! @param SP_CHAIN ShPtr to Chain to be inserted
    void ProteinModel::Insert( const util::ShPtr< Chain> &SP_CHAIN)
    {
      for
      (
        util::ShPtrVector< Chain>::iterator
          chain_itr( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        if( ( *chain_itr)->GetChainID() == SP_CHAIN->GetChainID())
        {
          for
          (
            storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( SP_CHAIN->GetData().Begin()),
              sse_itr_end( SP_CHAIN->GetData().End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            ( *chain_itr)->Insert( *sse_itr);
          }

          // emit change signal
          m_ChangeSignal.Emit( *this);
          return;
        }
        if( ( *chain_itr)->GetChainID() > SP_CHAIN->GetChainID())
        {
          m_Chains.InsertElement( chain_itr, SP_CHAIN);

          // emit change signal
          m_ChangeSignal.Emit( *this);
          return;
        }
      }

      m_Chains.PushBack( SP_CHAIN);

      // emit change signal
      m_ChangeSignal.Emit( *this);
    }

    //! @brief adds a structure element to the model
    //! @brief SP_SSE ShPtr to the SSE to be inserted
    //! @return whether insertion succeeded
    bool ProteinModel::Insert( const util::ShPtr< SSE> &SP_SSE)
    {
      // iterate over all chains, and find the right Chain mathich sse's chain id
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // if matching chain id and sse is not already contained, insert
        if( ( *chain_itr)->GetChainID() == SP_SSE->GetChainID() && !( *chain_itr)->DoesContain( *SP_SSE))
        {
          util::ShPtr< Chain> new_chain( chain_itr->HardCopy());

          if( new_chain->Insert( SP_SSE))
          {
            ( *chain_itr) = new_chain;

            // emit change signal
            m_ChangeSignal.Emit( *this);

            return true;
          }
          else
          {
            return false;
          }
        }
      }

      // end
      BCL_MessageVrb
      (
        "could not insert SSE " + SP_SSE->GetIdentification() + ", since no matching Chain was found or overlapping sse exists"
      );
      return false;
    }

    //! @brief replaces given SSE if it was in one of the chains
    //! @param SP_SSE SSE of interest
    //! @return whether replacement succeeded
    bool ProteinModel::Replace( const util::ShPtr< SSE> &SP_SSE)
    {
      //iterate over all chains, and find the right Chain matching sse's chain id
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // if matching chain id and if similar sse is within chain
        if( ( *chain_itr)->GetChainID() == SP_SSE->GetChainID() && ( *chain_itr)->DoesContain( *SP_SSE))
        {
          // replace chain and sse in chain
          util::ShPtr< Chain> sp_chain_new( chain_itr->HardCopy());

          if( sp_chain_new->Replace( SP_SSE))
          {
            ( *chain_itr) = sp_chain_new;

            // emit change signal
            m_ChangeSignal.Emit( *this);

            return true;
          }
        }
      }

      //end
      BCL_MessageCrt
      (
        "fail to replace SSELEMENT, " + SP_SSE->GetIdentification() +
        " since no matching Chain found in protein model with chains " + util::Format()( GetChainIDs())
      );
      return false;
    }

    //! @brief replaces given chain if it was in the protein model
    //! @param SP_CHAIN ShPtr to the chain of interest
    //! @return whether replacement succeeded
    bool ProteinModel::Replace( const util::ShPtr< Chain> &SP_CHAIN)
    {
      // iterate over all chains to find the right chain
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin());
          chain_itr != m_Chains.End();
        ++chain_itr
      )
      {
        // find the matching chain
        if( ( *chain_itr)->GetChainID() == SP_CHAIN->GetChainID())
        {
          // replace the chain
          ( *chain_itr) = SP_CHAIN;

          // emit change signal
          m_ChangeSignal.Emit( *this);

          // return true that the replacement was successful
          return true;
        }
      }

      // return false if replacement failed
      BCL_MessageCrt
      (
        "fail to replace CHAIN, " + SP_CHAIN->GetClassIdentifier() +
        " since no matching Chain found in protein model with chains " + util::Format()( GetChainIDs())
      );
      return false;
    }

    //! @brief replaces SSEs in the protein model with the given SSE
    //! @detail SSEs in the protein model that are overlapping with the given SSE are either completely
    //! replaced or resized in order to not overlap with the given SSE
    //! @param SP_SSE SSE that is to be inserted into the protein model
    void ProteinModel::ReplaceResize( const util::ShPtr< SSE> &SP_SSE)
    {
      // determine SSEs that are overlapping with the given SSE
      const util::SiPtrList< const SSE> overlap_sses( GetOverlappingSSEs( *SP_SSE));

      // determine the neighbor SSEs of the given SSE
      storage::VectorND< 2, util::SiPtr< const SSE> > neighbors( GetNeighborSSEs( *SP_SSE));
      if( !neighbors.First().IsDefined() || !neighbors.Second().IsDefined())
      {
        return;
      }

      const SSE neighbor_n( *neighbors.First());
      const SSE neighbor_c( *neighbors.Second());

      // determine if the neighbor SSEs have to be resized
      const int seq_n_n( neighbor_n.GetFirstAA()->GetSeqID());
      const int seq_c_n( neighbor_n.GetLastAA()->GetSeqID());
      const int seq_n_new( SP_SSE->GetFirstAA()->GetSeqID());
      const int seq_c_new( SP_SSE->GetLastAA()->GetSeqID());
      const int seq_n_c( neighbor_c.GetFirstAA()->GetSeqID());
      const size_t distance_n( seq_n_new - seq_c_n);
      const size_t distance_c( seq_n_c - seq_c_new);

      ReplaceWithOverlapping( SP_SSE);

      // resize n-terminal neighbor SSE if needed
      if( distance_n > 1)
      {
        const size_t length_n( neighbor_n.GetSize() + distance_n - 1);
        const biol::AASequence &sequence( *GetChain( neighbor_n.GetChainID())->GetSequence());
        const util::ShPtr< SSE> sp_new_n
        (
          new SSE( sequence.SubSequence( seq_n_n - 1, length_n), neighbor_c.GetType())
        );
        ReplaceWithOverlapping( sp_new_n);
      }
      else if( distance_n < 1)
      {
        BCL_MessageTop( "Shrinking n-terminal SSE.");
      }

      // resize c-terminal neighbor SSE if needed
      if( distance_c > 1)
      {
        const size_t length_c( neighbor_c.GetSize() + distance_c - 1);
        BCL_MessageTop
        (
          "Extending c-terminal SSE to " + util::Format()( seq_c_new + 1) + "-" +
          util::Format()( seq_c_new + length_c)
        );
        const biol::AASequence &sequence( *GetChain( neighbor_c.GetChainID())->GetSequence());
        const util::ShPtr< SSE> sp_new_c
        (
          new SSE( sequence.SubSequence( seq_c_new, length_c), neighbor_c.GetType())
        );
        ReplaceWithOverlapping( sp_new_c);
      }
      else if( distance_c < 1)
      {
        BCL_MessageTop( "Shrinking c-terminal SSE.");
      }
    }

    //! @brief replaces vector of SSEs in the chains
    //! @brief SSELEMENTS ShPtrVector of SSEs to be replaced
    //! @return whether replacement succeeded
    bool ProteinModel::Replace( const util::ShPtrVector< SSE> &SSELEMENTS)
    {
      //iterate over all given sses and replace them individually
      for
      (
        util::ShPtrVector< SSE>::const_iterator sse_itr( SSELEMENTS.Begin()),
          sse_itr_end( SSELEMENTS.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // try to replace single given sse
        if( !Replace( *sse_itr))
        {
          BCL_Exit( "unable to replace sselement", -1);
          return false;
        }
      }

      // no onchange signal emission necessary, since the Replace already emits

      //end
      return true;
    }

    //! @brief replace all SSEs that overlap with SP_SSE with SP_SSE
    //! @param SP_SSE ShPtr to SSE to be inserted
    //! @return whether replacement succeeded
    bool ProteinModel::ReplaceWithOverlapping( const util::ShPtr< SSE> &SP_SSE)
    {
      //iterate over all chains, and find the right Chain matching sse's chain id
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // if chain ids match and chain contains overlapping sses to SP_SSE, replace chain and sse in chain
        if
        (
          ( *chain_itr)->GetChainID() == SP_SSE->GetChainID() &&
          ( *chain_itr)->DoesContainOverlapping( *SP_SSE)
        )
        {
          util::ShPtr< Chain> sp_chain_new( ( *chain_itr)->Clone());
          ( *chain_itr) = sp_chain_new;

          // emit change signal
          m_ChangeSignal.Emit( *this);

          return sp_chain_new->ReplaceWithOverlapping( SP_SSE);
        }
      }

      // end
      BCL_MessageCrt
      (
        "fail to replace SSELEMENT " + SP_SSE->GetIdentification() +
        ", since no matching Chain or overlapping sse was found"
      );

      return false;
    }

    //! @brief replaces given set of SSEs in the chains
    //! @param SSELEMENTS set of SSEs to be replaced
    //! @return whether replacement succeeded
    bool ProteinModel::Replace( const storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> &SSELEMENTS)
    {
      //iterate over all given sses and replace them individually
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( SSELEMENTS.Begin()),
          sse_itr_end( SSELEMENTS.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        //try to replace single given sse
        if( !Replace( *sse_itr))
        {
          BCL_Exit( "unable to replace sselement", -1);
          return false;
        }
      }

      // no onchange signal emission necessary, since the Replace already emits

      //end
      return true;
    }

    //! @brief remove given SSE
    //! @brief THIS_SSE SSE of interest
    //! @return whether removal succeeded
    bool ProteinModel::Remove( const SSE &THIS_SSE)
    {
      //iterate over all chains, and find the right Chain matching sse's chain id
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // if matching chain id and if similar sse is within chain
        if( ( *chain_itr)->GetChainID() == THIS_SSE.GetChainID() && ( *chain_itr)->DoesContain( THIS_SSE))
        {
          util::ShPtr< Chain> new_chain( chain_itr->HardCopy());

          if( new_chain->Remove( THIS_SSE))
          {
            // replace chain
            ( *chain_itr) = new_chain;

            // emit change signal
            m_ChangeSignal.Emit( *this);

            return true;
          }
          else
          {
            return false;
          }
        }
      }

      //end
      BCL_MessageCrt
      (
        "fail to remove SSELEMENT: " + THIS_SSE.GetIdentification() + ", since no matching SSE found in protein model"
      );
      return false;
    }

    //! @brief sets all positions of all SSEs to ideal conformation and superimpose them with prior coordinates
    //! depending on flag KEEP_POSITION
    //! @param KEEP_POSITION flag to determine superimposition with prior coordinates
    void ProteinModel::SetToIdealConformation( const bool KEEP_POSITION)
    {
      //iterate over all chains
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        //set each chain to ideal conformation
        ( *chain_itr)->SetToIdealConformation( KEEP_POSITION);
      }

      // emit change signal
      m_ChangeSignal.Emit( *this);
    }

    //! @brief chop all sselements of that model in pieces of the given size vector MIN_SSE_LENGTHS
    //! @param MIN_SSE_LENGTHS minimal sizes for chopping for each SSType
    void ProteinModel::ChopSSEs( const storage::VectorND< 3, size_t> &MIN_SSE_LENGTHS)
    {
      //instantiate new vector to chopped elements
      util::ShPtrVector< SSE> new_model;

      //iterate over all sselements
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // chop sses for this chain
        ( *chain_itr)->ChopSSEs( MIN_SSE_LENGTHS);
      }

      // emit change signal
      m_ChangeSignal.Emit( *this);
    }

    //! @brief checks whether the given SSE is already in the model or is overlapping with one
    //! @param THIS_SSE SSE of interest
    //! @return whether the given SSE is already in the model or is overlapping with one
    bool ProteinModel::DoesContain( const SSE &THIS_SSE) const
    {
      //iterate over all sselements
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        if( ( *chain_itr)->GetChainID() == THIS_SSE.GetChainID() && ( *chain_itr)->DoesContain( THIS_SSE))
        {
          return true;
        };
      }

      return false;
    }

    //! @brief returns pairs of sses that have short loops (at most MAX_LOOP_LENGTH) between each other
    //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
    //! @return list of pairs of sses that have short loops (at most MAX_LOOP_LENGTH) between
    //! @return each other
    storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > > ProteinModel::GetSSEsWithShortLoops
    (
      const size_t MAX_LOOP_LENGTH
    ) const
    {
      // initialize list of pair of sses to be returned
      storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > > sse_list;

      // iterate over chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // append the list of sses with short loops for this chain
        sse_list.Append( ( *chain_itr)->GetSSEsWithShortLoops( MAX_LOOP_LENGTH));
      }
      // end
      return sse_list;
    }

    //! @brief returns SiPtrList of sses that have short loops ( at most MAX_LOOP_LENGTH) to provided TARGET_SSE
    //! @param TARGET_SSE SSE for which short loop connecting SSEs are being search
    //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
    //! @return SiPtrList of sses that have short loops to provided TARGET_SSE
    util::SiPtrList< const SSE> ProteinModel::GetSSEsWithShortLoops
    (
      const SSE &TARGET_SSE,
      const size_t MAX_LOOP_LENGTH
    ) const
    {
      // get the chain with the same chain ID as TARGET_SSE
      util::ShPtr< Chain> chain_ptr( GetChain( TARGET_SSE.GetChainID()));

      // if chain was not found
      if( chain_ptr.IsDefined())
      {
        // return this chain overlapping sses
        return chain_ptr->GetSSEsWithShortLoops( TARGET_SSE, MAX_LOOP_LENGTH);
      }

      // otherwise if chain was not found warn user
      BCL_MessageStd
      (
        "No sses with short loops to given list of sses are found in the protein model!"
      );

      // return empty list
      return util::SiPtrList< const SSE>();
    }

    //! @brief selects from provided SSE_LIST, sses that have short loops ( <=MAX_LOOP_LENGTH) to SSEs in ProteinModel
    //! @param SSE_LIST list of SSEs on which the selection will be done
    //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
    //! @return Subset of SSE_LIST that has short loops to SSEs in this ProteinModel
    util::SiPtrList< const SSE> ProteinModel::GetSSEsWithShortLoops
    (
      const util::SiPtrList< const SSE> &SSE_LIST,
      const size_t MAX_LOOP_LENGTH
    ) const
    {
      // initialize the list to be returned
      util::SiPtrList< const SSE> eligible_sses;

      // iterate over every SSE in this list
      for
      (
        util::SiPtrList< const SSE>::const_iterator sse_itr( SSE_LIST.Begin()), sse_itr_end( SSE_LIST.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // if the associated chain has at least 1 sses with a short loop to SSE behind sse_itr
        if( !( GetChain( ( *sse_itr)->GetChainID())->GetSSEsWithShortLoops( **sse_itr, MAX_LOOP_LENGTH)).IsEmpty())
        {
          // insert this sse into eligible_sses
          eligible_sses.PushBack( *sse_itr);
        }
      }

      // return
      return eligible_sses;
    };

    //! @brief returns the SSE before and SSE after given TARGET_SSE in this ProteinModel
    //! @param TARGET_SSE SSE of interest
    //! @return the SSE before and SSE after given TARGET_SSE in this ProteinModel
    storage::VectorND< 2, util::SiPtr< const SSE> > ProteinModel::GetNeighborSSEs
    (
      const SSE &TARGET_SSE
    ) const
    {
      // initialize the list to be returned
      storage::VectorND< 2, util::SiPtr< const SSE> > neighbor_sses;

      // get the chain with the same chain ID as TARGET_SSE
      util::ShPtr< Chain> chain_ptr( GetChain( TARGET_SSE.GetChainID()));

      // if chain is found
      if( chain_ptr.IsDefined())
      {
        // return the neighbor sses from this chain
        return chain_ptr->GetNeighborSSEs( TARGET_SSE);
      }

      // otherwise if chain was not found warn user
      BCL_MessageStd
      (
        "No neighbor SSEs are found for the given SSE, since ProteinModel does not have the specified chain!"
        + std::string( 1, TARGET_SSE.GetChainID())
      );

      // return
      return neighbor_sses;
    }

    //! @brief returns SSEs in this ProteinModel which overlap with given TARGET_SSE
    //! @brief TARGET_SSE SSE for which overlaps are going to be searched
    //! @return SSEs in this ProteinModel which overlap with given TARGET_SSE
    util::SiPtrList< const SSE> ProteinModel::GetOverlappingSSEs
    (
      const SSE &TARGET_SSE
    ) const
    {
      // return the associated chains overlapping sses with this one
      return GetChain( TARGET_SSE.GetChainID())->GetOverlappingSSEs( TARGET_SSE);
    }

    //! @brief returns SSEs in this ProteinModel which overlap with one or more SSEs in the given SSE_LIST
    //! @brief SSE_LIST SiPtrList of SSEs for which overlaps are going to be searched
    //! @return SSEs in this ProteinModel which overlap with one or more SSEs in the given SSE_LIST
    util::SiPtrList< const SSE> ProteinModel::GetOverlappingSSEs
    (
      const util::SiPtrList< const SSE> &SSE_LIST
    ) const
    {
      // initialize the list to be returned
      util::SiPtrList< const SSE> eligible_sses;

      // iterate over every SSE in this list
      for
      (
        util::SiPtrList< const SSE>::const_iterator sse_itr( SSE_LIST.Begin()), sse_itr_end( SSE_LIST.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // get overlapping sses for this sse
        util::SiPtrList< const SSE> overlapping_sses( GetOverlappingSSEs( **sse_itr));

        // if any were found
        if( !overlapping_sses.IsEmpty())
        {
          // pushback overlapping sses in the model for this sse into the list
          eligible_sses.Append( GetOverlappingSSEs( **sse_itr));
        }
      }

      // return
      return eligible_sses;
    }

    //! @brief AddLoops generates loop SSEs for the protein model and assigns them zero coordinates
    //! @param UNDEFINED_COORDINATES create loop with undefined coordinates, or with coordinates form the member sequence
    //! @param MERGE_CONSECUTIVE_SSES merge consecutive sses of given type
    //! @param SS_TYPE hte sstype of the sses to be merged
    void ProteinModel::AddLoops
    (
      const bool UNDEFINED_COORDINATES,
      const bool MERGE_CONSECUTIVE_SSES,
      const biol::SSType &SS_TYPE
    )
    {
      // iterate through all chains
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        ( *chain_itr)->AddLoops( UNDEFINED_COORDINATES, MERGE_CONSECUTIVE_SSES, SS_TYPE);
      }
    }

    //! @brief filters the current protein model by given minimum SSE sizes
    //! @param MIN_SSE_SIZES minimum SSE sizes to filter the model by
    void ProteinModel::FilterByMinSSESizes
    (
      const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES
    )
    {
      // iterate through all chains
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // make a clone of the chain and call the function
        ( *chain_itr) = util::ShPtr< Chain>( ( *chain_itr)->Clone());
        ( *chain_itr)->FilterByMinSSESizes( MIN_SSE_SIZES);
      }

      // emit change signal
      m_ChangeSignal.Emit( *this);
    }

    //! @brief join following ( progressing sequence id) SSEs of given SS_TYPE into one SSE
    //! @param SS_TYPE SSType of interest
    //! @param TEST_PEPTIDE_BOND join only, if SSEs are connected by peptide bond
    void ProteinModel::Join( const biol::SSType &SS_TYPE, const bool TEST_PEPTIDE_BOND)
    {
      // iterate through all chains
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // make a copy of the chain and call the function
        ( *chain_itr) = util::ShPtr< Chain>( ( *chain_itr)->HardCopy());
        ( *chain_itr)->Join( SS_TYPE, TEST_PEPTIDE_BOND);
      }

      // emit change signal
      m_ChangeSignal.Emit( *this);
    }

    //! @brief let the aadata of the each sse point to the corresponding data in the chain - determined by pdbID
    void ProteinModel::ConnectSSEToChainData()
    {
      // iterate through the chains
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        ( *chain_itr)->ConnectSSEToChainData();
      }

      // emit change signal
      m_ChangeSignal.Emit( *this);
    }

    //! @brief Replace SSEs with those drawn from the pool
    void ProteinModel::AdoptSSEsMaintainCoordinates( const util::SiPtrVector< const SSE> &SSES)
    {
      storage::Map< char, util::SiPtrVector< const SSE> > chain_to_sses;
      for( auto itr_sses( SSES.Begin()), itr_sses_end( SSES.End()); itr_sses != itr_sses_end; ++itr_sses)
      {
        chain_to_sses[ ( *itr_sses)->GetChainID()].PushBack( *itr_sses);
      }

      // iterate through the chains
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        ( *chain_itr)->AdoptSSEsMaintainCoordinates( chain_to_sses[ ( *chain_itr)->GetChainID()]);
      }

      // emit change signal
      m_ChangeSignal.Emit( *this);
    }

    //! @brief Get SSE hash string to aid in identifying similar proteins
    std::string ProteinModel::GetSSEHashString() const
    {
      std::string hash;
      // iterate through the chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        hash += ( *chain_itr)->GetChainID();
        hash += '{';
        hash += ( *chain_itr)->GetSSEHashString();
        hash += '}';
      }
      return hash;
    }

    //! @brief calculates all pairwise interactions of amino acids within SSE fragments
    //! @param SIZE number of amino acids in the chain
    //! @return matrix of pairwise interactions of amino acids within SSE fragments
    const linal::Matrix< double> ProteinModel::CalculateFragmentInteractions() const
    {
      // get size of matrix
      const size_t size( GetNumberAAs());

      // initialize matrix
      linal::Matrix< double> matrix( size, size);

      // get fragments
      util::ShPtrVector< SSEGeometryInterface> fragments( GetFragments());

      // criteria to determine whether or not two fragments are in contact
      const math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> &PACKING_CRITERIA( *CollectorTopologyCombined::GetPackingCriteria());

      //get ids of fragments
      storage::Vector< storage::Pair< int, int> > fragment_ids( GetFragmentIDs());

      // get sses
      const util::SiPtrVector< const SSE> sses( GetSSEs());

      // initialize weight of fragment interaction
      double interaction_weight( 0);

      //iterate over ids to fill in row of table

      // initialize fragment vector for each sse
      util::ShPtrVector< SSEGeometryInterface> sse_fragments;

      for
      (
        util::ShPtrVector< SSEGeometryInterface>::const_iterator frag_itr( fragments.Begin()),
        frag_itr_end( fragments.End()); frag_itr != frag_itr_end; ++frag_itr
      )
      {
        for
        (
          util::ShPtrVector< SSEGeometryInterface>::const_iterator frag_itr_2( fragments.Begin() + 1),
          frag_itr_2_end( fragments.End()); frag_itr_2 != frag_itr_2_end; ++frag_itr_2
        )
        {
          //get packing between the two fagments
          SSEGeometryPacking this_packing( ( *GetSSEGeometryPackingPickers().e_BestInteractionWeight)->operator()( **frag_itr, **frag_itr_2));
          if( PACKING_CRITERIA( this_packing))
          {
            //fill in matrix with interaction weight as metric of packing score
            interaction_weight = this_packing.GetInteractionWeight();
            matrix( ( *frag_itr_2)->GetCentralAA(), ( *frag_itr)->GetCentralAA()) = interaction_weight;
          }
        }
      }

      return matrix;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    //! @param PROTEIN_MODEL ProteinModel to be copied
    //! @return This model after all members are assigned to values from PROTEIN_MODEL
    ProteinModel &ProteinModel::operator =( const ProteinModel &PROTEIN_MODEL)
    {
      // update members
      m_Chains = PROTEIN_MODEL.m_Chains;
      m_ProteinModelData = PROTEIN_MODEL.m_ProteinModelData;

      // emit change signal
      m_ChangeSignal.Emit( *this);

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read ProteinModel from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModel::Read( std::istream &ISTREAM)
    {
      //read data
      io::Serialize::Read( m_Chains, ISTREAM);
      io::Serialize::Read( m_ProteinModelData, ISTREAM);

      // emit change signal
      m_ChangeSignal.Emit( *this);

      //return
      return ISTREAM;
    }

    //! @brief write ProteinModel to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &ProteinModel::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      //write data
      io::Serialize::Write( m_Chains, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ProteinModelData, OSTREAM, INDENT);

      //return
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
