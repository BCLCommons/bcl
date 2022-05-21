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
#include "contact/bcl_contact_map.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

  //////////
  // data //
  //////////

    //! identifier string to be used for prediction maps
    const std::string Map::s_Identifier = "CHAINS";

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Map::Map() :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< size_t>, bool> >(),
      m_Chains(),
      m_Boundary( 8)
    {
    }

    //! @brief construct contact map from a chain
    //! @param CHAIN Chain for which sequence contacts are going to be deteced
    //! @param BOUNDARY number of residues to be excluded from both sides
    Map::Map
    (
      const util::ShPtr< assemble::Chain> &CHAIN,
      const int BOUNDARY
    ) :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< size_t>, bool> >(),
      m_Chains( 1, CHAIN),
      m_Boundary( BOUNDARY)
    {
      FillMap();
    }

    //! @brief construct contact map from a protein model( intra-sequence contacts and inter-sequence contacts),
    //! @param PROTEIN_MODEL Protein Model for which inter and intra sequence contacts are going to be deteced
    //! @param BOUNDARY number of residues to be excluded from both sides
    Map::Map
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const int BOUNDARY
    ) :
      storage::ObjectNDHashMap< 2, biol::AAData, storage::Pair< linal::Vector< size_t>, bool> >(),
      m_Chains( PROTEIN_MODEL.GetChains()),
      m_Boundary( BOUNDARY)
    {
      FillMap();
    }

    //! @brief virtual copy constructor
    Map *Map::Clone() const
    {
      return new Map( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Map::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns stored chains
    //! @return stored chains
    const util::ShPtrVector< assemble::Chain> &Map::GetChains() const
    {
      return m_Chains;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief iterates over the stored chains and finds the one with matching CHAIN_ID
    //! @param CHAIN_ID chain id of the chain that is beings searched
    //! @return Chain  with the searched chain id
    const util::ShPtr< assemble::Chain> &Map::GetChain( const char CHAIN_ID) const
    {
      // iterate over every chain stored
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // if the sequence id matches return it
        if( ( *chain_itr)->GetSequence()->GetChainID() == CHAIN_ID)
        {
          return *chain_itr;
        }
      }

      // else Exit
      BCL_Exit( "No chain with the provided chain id \'" + util::Format()( CHAIN_ID) + "\' is not stored!!", -1);
      return m_Chains.FirstElement();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Returns contacts vector for the provided AA Data pointers
    //! @param AA_DATA_POINTERS pair of pointers to AAData of amino acids of interest
    //! @return the pair of vector of contact vectors and the merged prediction
    const storage::Pair< linal::Vector< size_t>, bool> &Map::GetContactVector
    (
      const storage::VectorND< 2, util::SiPtr< const biol::AAData> > &AA_DATA_POINTERS
    ) const
    {
      // initialize undefined predictions vector
      static const storage::Pair< linal::Vector< size_t>, bool> s_undefined_contact_vector
      (
        linal::Vector< size_t>( Types::s_NumberValidTypes, util::GetUndefined< size_t>()),
        false
      );

      // search in the hash map
      storage::HashMap< size_t, storage::Pair< linal::Vector< size_t>, bool> >::const_iterator itr
      (
        Find( AA_DATA_POINTERS)
      );

      // if not found return undefined
      if( itr == End())
      {
        return s_undefined_contact_vector;
      }

      // end
      return itr->second;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Map::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &Map::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    //! @brief helper function to read predictions from a file
    //! @param ISTREAM input stream to be read from
    //! @return input stream from which the map was read
    std::istream &Map::ReadMap( std::istream &ISTREAM)
    {
      // header that is going to be searched
      std::string header;

      // while reading header and not end of file yet
      while( ISTREAM >> header && !ISTREAM.eof())
      {
        // if header matches the identifier
        if( header == s_Identifier)
        {

          // read the two chain ids
          char chain_id_a, chain_id_b;
          ISTREAM >> chain_id_a >> chain_id_b;

          // check the chain_id_a exists
          BCL_Assert
          (
            GetChain( chain_id_a).IsDefined(),
            "No chain is stored in the map with given chain id: " + std::string( 1, chain_id_a)
          );

          // check the chain_id_b exists
          BCL_Assert
          (
            GetChain( chain_id_b).IsDefined(),
            "No chain is stored in the map with given chain id: " + std::string( 1, chain_id_b)
          );

          // read the contact map for two chains with the provided chain ids
          ReadContacts
          (
            ISTREAM,
            *GetChain( chain_id_a),
            *GetChain( chain_id_b)
          );
        }
      }

      // end
      return ISTREAM;
    }

    //! @brief helper function to write predictions to a file
    //! @param OSTREAM output stream to write to
    //! @param WRITE_ONLY_CONTACTS boolean value that determines whether non contacting residue couples are written
    //! @return output stream which was written to
    std::ostream &Map::WriteMap
    (
      std::ostream &OSTREAM,
      const bool WRITE_ONLY_CONTACTS
    ) const
    {
      // iterate over every sequence
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr_a( m_Chains.Begin()),
        chain_itr_end( m_Chains.End());
        chain_itr_a != chain_itr_end;
        ++chain_itr_a
      )
      {
        // versus every other sequence
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr_b( m_Chains.Begin());
          chain_itr_b != chain_itr_end;
          ++chain_itr_b
        )
        {
          // output the header
          OSTREAM << s_Identifier << ' '
                  << util::Format()( (*chain_itr_a)->GetChainID()) << ' '
                  << util::Format()( (*chain_itr_b)->GetChainID()) << '\n';

          // write predictions
          WriteContacts( OSTREAM, **chain_itr_a, **chain_itr_b, WRITE_ONLY_CONTACTS);
        }
      }
      return OSTREAM;
    }

    //! @brief reads contacts for the provided chains
    //! @param ISTREAM input stream
    //! @param CHAIN_A first chain
    //! @param CHAIN_B second chain
    //! @return istream from which was read
    std::istream &Map::ReadContacts
    (
      std::istream &ISTREAM,
      const assemble::Chain &CHAIN_A,
      const assemble::Chain &CHAIN_B
    )
    {
      // initialize sequence id and AAType pairs and the vector hold predictions
      int seq_id_a, seq_id_b;
      char type_a, type_b;
      storage::Pair< linal::Vector< size_t>, bool> contact_vector
      (
        linal::Vector< size_t>( Types::s_NumberValidTypes, util::GetUndefined< size_t>()),
        false
      );

      // while it is possible to read
      while( ISTREAM >> seq_id_a >> type_a >> seq_id_b >> type_b && !ISTREAM.eof())
      {
        // make the sure sequences match
        BCL_Assert
        (
          CHAIN_A.GetSequence()->GetAA( seq_id_a - 1)->GetType()->GetOneLetterCode() == type_a &&
          CHAIN_B.GetSequence()->GetAA( seq_id_b - 1)->GetType()->GetOneLetterCode() == type_b,
          " The provided contact pair does not match the provided sequences" +
          util::Format()( CHAIN_A.GetSequence()->GetAA( seq_id_a - 1)->GetType()->GetOneLetterCode()) +
          " vs " + util::Format()( type_a) + " and " +
          util::Format()( CHAIN_B.GetSequence()->GetAA( seq_id_b - 1)->GetType()->GetOneLetterCode()) +
          " vs " + util::Format()( type_b)
        );

        // read the predictions vector
        ISTREAM >> contact_vector.First()( GetTypes().HELIX_HELIX)
                >> contact_vector.First()( GetTypes().HELIX_SHEET)
                >> contact_vector.First()( GetTypes().SHEET_HELIX)
                >> contact_vector.First()( GetTypes().STRAND_STRAND)
                >> contact_vector.First()( GetTypes().SHEET_SHEET);

        // read the is_in_contact from the file
        ISTREAM >> contact_vector.Second();

        // since default constructed conteact vectors have 0, but when reading all values are set to undefined
        // HELIX_STRAND and STRAND_HELIX contacts cannot be read and will remain undefined
        contact_vector.First()( GetTypes().HELIX_STRAND) = 0;
        contact_vector.First()( GetTypes().STRAND_HELIX) = 0;

        // Insert amino acid pair
        Insert
        (
          storage::VectorND< 2, util::SiPtr< const biol::AAData> >
          (
            CHAIN_A.GetSequence()->GetAA( seq_id_a - 1)->GetData(),
            CHAIN_B.GetSequence()->GetAA( seq_id_b - 1)->GetData()
          ),
          contact_vector
        );
      }

      // end
      return ISTREAM;
    }

    //! @brief Writes predictions for the provided chains
    //! @param OSTREAM output stream
    //! @param CHAIN_A first chain
    //! @param CHAIN_B second chain
    //! @param WRITE_ONLY_CONTACTS boolean value that determines whether non contacting residue couples are written
    //! @return ostream which was written to
    std::ostream &Map::WriteContacts
    (
      std::ostream &OSTREAM,
      const assemble::Chain &CHAIN_A,
      const assemble::Chain &CHAIN_B,
      const bool WRITE_ONLY_CONTACTS
    ) const
    {
      // iterate over first sequence
      for
      (
        biol::AASequence::const_iterator aa_itr_a( CHAIN_A.GetSequence()->Begin()),
          aa_itr_a_end( CHAIN_A.GetSequence()->End());
        aa_itr_a != aa_itr_a_end;
        ++aa_itr_a
      )
      {
        // iterate over second sequence
        for
        (
          biol::AASequence::const_iterator aa_itr_b( CHAIN_B.GetSequence()->Begin()),
            aa_itr_b_end( CHAIN_B.GetSequence()->End());
          aa_itr_b != aa_itr_b_end;
          ++aa_itr_b
        )
        {
          // get the contact vector for the two amino acids
          storage::Pair< linal::Vector< size_t>, bool> contact_vector
          (
            GetContactVector
            (
              storage::VectorND< 2, util::SiPtr< const biol::AAData> >
              (
                ( *aa_itr_a)->GetData(),
                ( *aa_itr_b)->GetData()
              )
            )
          );

          // if these residues are not in contact and WRITE_ONLY_CONTACTS flag is set skip
          if( WRITE_ONLY_CONTACTS && contact_vector.Second() == 0)
          {
            continue;
          }

          // if these residues are not found
          if( !contact_vector.First().IsDefined())
          {
            // if write_only_contacts flag is true skip
            if( WRITE_ONLY_CONTACTS)
            {
              continue;
            }
            else
            {
              contact_vector.First() = linal::Vector< size_t>( size_t( Types::s_NumberValidTypes), size_t( 0));
              contact_vector.Second() = false;
            }
          }

          // output values
          OSTREAM << util::Format().W( 4)( ( *aa_itr_a)->GetSeqID()) << ' '
                  << ( *aa_itr_a)->GetType()->GetOneLetterCode() << ' '
                  << util::Format().W( 4)( ( *aa_itr_b)->GetSeqID()) << ' '
                  << ( *aa_itr_b)->GetType()->GetOneLetterCode() << ' '
                  << util::Format()( contact_vector.First()( GetTypes().HELIX_HELIX)) << ' '
                  << util::Format()( contact_vector.First()( GetTypes().HELIX_SHEET)) << ' '
                  << util::Format()( contact_vector.First()( GetTypes().SHEET_HELIX)) << ' '
                  << util::Format()( contact_vector.First()( GetTypes().STRAND_STRAND)) << ' '
                  << util::Format()( contact_vector.First()( GetTypes().SHEET_SHEET)) << ' '
                  << util::Format()( contact_vector.Second()) << '\n';
        }
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Given 2 AAs, this function returns a bool of whether these two residues
    //! @brief are in contact according to rules of the specified contact type
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return boolean indicating whether AMINO_ACID_A and AMINO_ACID_B are in contact with any contact type
    bool Map::IsInContact( const biol::AABase &AMINO_ACID_A, const biol::AABase &AMINO_ACID_B)
    {
      for( Types::const_iterator itr( GetTypes().Begin()), itr_end( GetTypes().End()); itr != itr_end; ++itr)
      {
        bool in_contact( IsInContact( AMINO_ACID_A, AMINO_ACID_B, *itr));
        if( in_contact)
        {
          return true;
        }
      }
      return false;
    }

    //! @brief Given 2 AAs and a contact type, this function returns a bool of whether these two residues
    //! @brief are in contact according to rules of the specified contact type
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @param CONTACT_TYPE contact type
    //! @return boolean indicating whether AMINO_ACID_A and AMINO_ACID_B are in contact with CONTACT_TYPE
    bool Map::IsInContact
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const Type &CONTACT_TYPE
    )
    {
      // if it is one of the valid contacttypes
      if( CONTACT_TYPE->IsValid())
      {
        // calculate the distance between CB atoms of two residues
        const double distance( Distance( AMINO_ACID_A.GetFirstSidechainAtom(), AMINO_ACID_B.GetFirstSidechainAtom()));

        // if  distance is smaller than the supplied threshold for that given contacttype
        if( distance < CONTACT_TYPE->GetResidueDistanceCutoff())
        {
          return true;
        }
      }
      // flag is supplied, check the unknown contacts
      if( CONTACT_TYPE == GetTypes().e_Undefined)
      {
        // calculate the distance between CB atoms of two residues
        const double distance( Distance( AMINO_ACID_A.GetFirstSidechainAtom(), AMINO_ACID_B.GetFirstSidechainAtom()));

        // if distance is smaller than the supplied threshold for undefined contacts
        if( distance < Types::GetUnknownResidueDistanceCutoff())
        {
          return true;
        }
      }
      // else return false
      return false;
    }

    //! @brief checks if the two SSEs are in contact in the given contact map
    //! @param CONTACTS the contact map the defines the contacts between amino acids
    //! @param FIRST first SSE
    //! @param SECOND second SSE
    //! @return true if the SSEs are in contact
    bool Map::IsInContact( const Map &CONTACTS, const assemble::SSE &FIRST, const assemble::SSE &SECOND)
    {
      util::SiPtrVector< const biol::AABase> first_aas( FIRST.GetData()), second_aas( SECOND.GetData());
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          itr_aa_a( first_aas.Begin()), itr_aa_a_end( first_aas.End());
        itr_aa_a != itr_aa_a_end; ++itr_aa_a
      )
      {
        for
        (
          util::SiPtrVector< const biol::AABase>::const_iterator
            itr_aa_b( second_aas.Begin()), itr_aa_b_end( second_aas.End());
          itr_aa_b != itr_aa_b_end; ++itr_aa_b
        )
        {
          if( CONTACTS.IsInContact( **itr_aa_a, **itr_aa_b))
          {
            return true;
          }
        }
      }
      return false;
    }

    //! @brief This function fills m_Map with real contact information ( 0|1) for each residue couple from
    //! @brief each contacttype for the given set of sses
    void Map::FillMap()
    {
      BCL_MessageDbg( "Initializing the contactmap");

      // make sure all chains have SSE information
      CheckChainsForSSEs();

      // iterate over every chain stored
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr_a( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr_a != chain_itr_end;
        ++chain_itr_a
      )
      {
        FillMap( **chain_itr_a);

        // iterate over every other chain
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr_b( chain_itr_a + 1);
          chain_itr_b != chain_itr_end;
          ++chain_itr_b
        )
        {
          FillMap( **chain_itr_a, **chain_itr_b);
        }
      }
    }

    //! @brief This function fills m_Map with real contact information ( 0|1) for each residue couple from
    //! @brief each contact type for the given chain
    //! @param CHAIN a chain
    void Map::FillMap( const assemble::Chain &CHAIN)
    {
      // iterate over every sse element( sse_itr_a)
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr_a( CHAIN.GetData().Begin()), sse_itr_end( CHAIN.GetData().End());
        sse_itr_a != sse_itr_end;
        ++sse_itr_a
      )
      {
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
        sse_itr_b( sse_itr_a);
        ++sse_itr_b;

        // iterate over every later sse element( sse_itr_b)
        for
        (
          ;
          sse_itr_b != sse_itr_end;
          ++sse_itr_b
        )
        {
          // precalculate the contact type between sses, since this will stay same for all aa in these two sses
          Type contact_type
          (
            assemble::SSEGeometryPacking( **sse_itr_a, **sse_itr_b).GetContactType()
          );

          // HELIX_STRAND and STRAND_HELIX not handled yet - set to undefined
          if( contact_type == GetTypes().HELIX_STRAND)
          {
            contact_type = GetTypes().UNDEFINED_HELIX_STRAND;
          }
          else if( contact_type == GetTypes().STRAND_HELIX)
          {
            contact_type = GetTypes().UNDEFINED_STRAND_HELIX;
          }

          BCL_MessageDbg( "Determined contact type :" + util::Format()( contact_type));

          // iterate over every residue in sse_itr_a
          for
          (
            biol::AASequence::const_iterator aa_itr_a( ( *sse_itr_a)->Begin()),
              aa_itr_end_a( ( *sse_itr_a)->End());
            aa_itr_a != aa_itr_end_a;
            ++aa_itr_a
          )
          {
            // vs every residue in sse_itr_b
            for
            (
              biol::AASequence::const_iterator aa_itr_b( ( *sse_itr_b)->Begin()),
                aa_itr_end_b( ( *sse_itr_b)->End());
              aa_itr_b != aa_itr_end_b;
              ++aa_itr_b
            )
            {
              // if lastseqid was set
              if( CheckBoundaryCondition( **aa_itr_a, CHAIN) && CheckBoundaryCondition( **aa_itr_b, CHAIN))
              {
                InsertAminoAcidPair( **aa_itr_a, **aa_itr_b, contact_type, true);
              }
            }
          }
        }
      }
    }

    //! @brief This function fills m_Map with real contact information ( 0|1) for each residue couple from
    //! @brief each contact type for the given pair of chains
    //! @param CHAIN_A first chain
    //! @param CHAIN_B second chain
    void Map::FillMap( const assemble::Chain &CHAIN_A, const assemble::Chain &CHAIN_B)
    {
      // iterate over every sse element( sse_itr_a) in first chain
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr_a( CHAIN_A.GetData().Begin()), sse_itr_end_a( CHAIN_A.GetData().End());
        sse_itr_a != sse_itr_end_a;
        ++sse_itr_a
      )
      {
        // iterate over every later sse element( sse_itr_b) in the second chain
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr_b( CHAIN_B.GetData().Begin()), sse_itr_end_b( CHAIN_B.GetData().End());
          sse_itr_b != sse_itr_end_b;
          ++sse_itr_b
        )
        {
          // precalculate the contactype between sses, since this will stay same for all aa in these two sses
          Type contact_type
          (
            assemble::SSEGeometryPacking( **sse_itr_a, **sse_itr_b).GetContactType()
          );

          // HELIX_STRAND and STRAND_HELIX not handled yet - set to undefined
          if( contact_type == GetTypes().HELIX_STRAND || contact_type == GetTypes().STRAND_HELIX)
          {
            contact_type = GetTypes().e_Undefined;
          }

          BCL_MessageDbg( "Determined contact type :" + util::Format()( contact_type));

          // iterate over every residue in sse_itr_a
          for
          (
            biol::AASequence::const_iterator aa_itr_a( ( *sse_itr_a)->Begin()),
              aa_itr_end_a( ( *sse_itr_a)->End());
            aa_itr_a != aa_itr_end_a;
            ++aa_itr_a
          )
          {
            // vs every residue in sse_itr_b
            for
            (
              biol::AASequence::const_iterator aa_itr_b( ( *sse_itr_b)->Begin()),
                aa_itr_end_b( ( *sse_itr_b)->End());
              aa_itr_b != aa_itr_end_b;
              ++aa_itr_b
            )
            {
              // if lastseqid was set
              if( CheckBoundaryCondition( **aa_itr_a, CHAIN_A) && CheckBoundaryCondition( **aa_itr_b, CHAIN_B))
              {
                InsertAminoAcidPair( **aa_itr_a, **aa_itr_b, contact_type, true);
              }
            }
          }
        }
      }
    }

    //! @brief checks the boundary conditions of a residue to determine whether it should be included in the map
    //! @param AMINO_ACID amino acid for which the conditions are going to be checked
    //! @param CHAIN
    bool Map::CheckBoundaryCondition( const biol::AABase &AMINO_ACID, const assemble::Chain &CHAIN) const
    {
      return
      (
        AMINO_ACID.GetSeqID() > m_Boundary
      )
      &&
      (
        AMINO_ACID.GetSeqID() < CHAIN.GetSequence()->GetLastAA()->GetSeqID() - m_Boundary
      );
    }

    //! @brief insert the provided amino acid pair and contact type and optional reverse contact type,
    //! @brief inserts the related data into hashmap and residue pair list
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @param CONTACT_TYPE contact type
    //! @param INSERT_REVERSE also insert reverse of CONTACT_TYPE
    void Map::InsertAminoAcidPair
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const Type &CONTACT_TYPE,
      const bool INSERT_REVERSE
    )
    {
      //! create the vectors for storage
      storage::Pair< linal::Vector< size_t>, bool> contact_vector
      (
        linal::Vector< size_t>( Types::s_NumberValidTypes, size_t( 0)),
        false
      );

      // determine if these two residues are in contact
      const bool is_in_contact( IsInContact( AMINO_ACID_A, AMINO_ACID_B, CONTACT_TYPE));

      // if in contact, update the vector
      if( is_in_contact && CONTACT_TYPE->IsValid())
      {
        contact_vector.First()( CONTACT_TYPE) = true;
      }

      // set the boolean value accordingly
      contact_vector.Second() = is_in_contact;

      // create simple pointers to amino acid data
      util::SiPtr< const biol::AAData> sp_a( AMINO_ACID_A.GetData());
      util::SiPtr< const biol::AAData> sp_b( AMINO_ACID_B.GetData());

      // insert h,i
      Insert
      (
        storage::VectorND< 2, util::SiPtr< const biol::AAData> >( sp_a, sp_b),
        contact_vector
      );

      // if reverse flag is set
      if( INSERT_REVERSE)
      {
        storage::Pair< linal::Vector< size_t>, bool> contact_vector_swapped( contact_vector);
        linal::Vector< size_t> &contact_vector( contact_vector_swapped.First());
        std::swap( contact_vector( GetTypes().HELIX_SHEET), contact_vector( GetTypes().SHEET_HELIX));

        // insert i,h
        Insert
        (
          storage::VectorND< 2, util::SiPtr< const biol::AAData> >( sp_b, sp_a),
          contact_vector_swapped
        );
      }
    }

    //! @brief this functions checks that the chains stored also have SSE information
    //! @brief so that the contact map can be generated correctly
    void Map::CheckChainsForSSEs() const
    {
      // iterate over every chain stored
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( m_Chains.Begin()),
          chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        //check that there is at least 1 sse stored with this chain
         BCL_Assert
         (
           !( *chain_itr)->GetData().IsEmpty(),
           "The following chains does not have any SSE information: " + util::Format()( ( *chain_itr)->GetChainID())
         );
      }
    }

  } // namespace contact
} // namespace bcl
