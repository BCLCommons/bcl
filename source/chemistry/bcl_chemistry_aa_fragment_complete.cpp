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
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "graph/bcl_graph_subgraph.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief create a map from biol::AAType, biol::AtomType to id of the atom type within the sdf::BondInfo
    //! @param C_TERMINAL true to create a map assuming the residue is the C-terminus
    //! @param N_TERMINAL true to create a map assuming the resiude is the N-terminus
    storage::Vector< storage::Vector< size_t> >
      CreateAATypeAtomTypeIndexer( const bool &C_TERMINAL, const bool &N_TERMINAL)
    {
      storage::Vector< storage::Vector< size_t> > map
      (
        biol::GetAATypes().GetEnumCount(),
        storage::Vector< size_t>( biol::GetAtomTypes().GetEnumCount(), util::GetUndefined< size_t>())
      );
      for
      (
        biol::AATypes::const_iterator itr_aa( biol::GetAATypes().Begin()), itr_aa_end( biol::GetAATypes().End());
        itr_aa != itr_aa_end;
        ++itr_aa
      )
      {
        const FragmentConfigurationShared &fragment( ( *itr_aa)->GetFragment( C_TERMINAL, N_TERMINAL));
        // get the property that indicates the biol::AtomTypes
        storage::Vector< std::string> types( util::SplitString( fragment.GetMDLProperty( "BiolAtomTypes"), " \n"));
        for( size_t atom_index( 0), max_atom_index( types.GetSize()); atom_index < max_atom_index; ++atom_index)
        {
          biol::AtomType atom_type( util::TrimString( types( atom_index)));
          if( atom_type.IsDefined())
          {
            map( itr_aa->GetIndex())( atom_type.GetIndex()) = atom_index;
          }
        }
      }
      return map;
    }

    //! @brief get a map from biol::AAType, biol::AtomType to id of the atom type within the sdf::BondInfo
    //! @param C_TERMINAL true to create a map assuming the residue is the C-terminus
    //! @param N_TERMINAL true to create a map assuming the resiude is the N-terminus
    const storage::Vector< storage::Vector< size_t> > &
      GetAATypeAtomTypeToIndexMap( const bool &C_TERMINAL, const bool &N_TERMINAL)
    {
      static storage::Vector< storage::Vector< size_t> >
        s_NCTerminal( CreateAATypeAtomTypeIndexer( true, true)),
        s_NTerminal( CreateAATypeAtomTypeIndexer( false, true)),
        s_CTerminal( CreateAATypeAtomTypeIndexer( true, false)),
        s_NonTerminal( CreateAATypeAtomTypeIndexer( false, false));
      return C_TERMINAL ? ( N_TERMINAL ? s_NCTerminal : s_CTerminal) : ( N_TERMINAL ? s_NTerminal : s_NonTerminal);
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAFragmentComplete::s_Instance
    (
      GetObjectInstances().AddInstance( new AAFragmentComplete())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    AAFragmentComplete::AAFragmentComplete()
    {
    }

    //! @brief construct from an AA sequence.  Normally the contained AA class should be AA complete
    //! @param ALLOW_UNDEFINED_POSITIONS true to allow undefined atom positions
    //! Default is to assert and fail if heavy atom position is absent
    AAFragmentComplete::AAFragmentComplete( const biol::AASequence &SEQUENCE, const bool &ALLOW_UNDEFINED_POSITIONS)
    {
      Construct( util::SiPtrVector< const biol::AABase>( SEQUENCE.Begin(), SEQUENCE.End()), ALLOW_UNDEFINED_POSITIONS);
    }

    //! @brief construct from an AA sequence.  Normally the contained AA type should be AA complete
    //! @param ALLOW_UNDEFINED_POSITIONS true to allow undefined atom positions
    //! Default is to assert and fail if heavy atom position is absent
    AAFragmentComplete::AAFragmentComplete
    (
      const util::SiPtrVector< const biol::AABase> &AA_SEQUENCE,
      const bool &ALLOW_UNDEFINED_POSITIONS
    )
    {
      Construct( AA_SEQUENCE, ALLOW_UNDEFINED_POSITIONS);
    }

    //! @brief Clone function
    //! @return pointer to new AAFragmentComplete
    AAFragmentComplete *AAFragmentComplete::Clone() const
    {
      return new AAFragmentComplete( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAFragmentComplete::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return each residue as a separate fragment
    //! @param CONSIDER_BACK_BONE whether to include the back bone atoms
    //! @param CONSIDER_SIDE_CHAIN whether to include the side chain atoms
    FragmentEnsemble AAFragmentComplete::GetResiduesAsFragments
    (
      const bool &CONSIDER_BACK_BONE,
      const bool &CONSIDER_SIDE_CHAIN
    ) const
    {
      BCL_Assert
      (
        CONSIDER_BACK_BONE || CONSIDER_SIDE_CHAIN,
        "Must consider either or both the back bone and/or side chain"
      );

      FragmentEnsemble ensemble;

      // indices of all atoms belonging to each residue
      storage::Vector< storage::Vector< size_t> > residue_atom_indices( m_Sequence.GetSize());

      // determine indices of atoms for each residue
      size_t atom_index( 0);
      for
      (
        storage::Vector< size_t>::const_iterator itr( m_AAIndices.Begin()), itr_end( m_AAIndices.End());
        itr != itr_end;
        ++itr, ++atom_index
      )
      {
        if( m_BiolAtomType( atom_index)->IsSideChain() == CONSIDER_SIDE_CHAIN || CONSIDER_BACK_BONE)
        {
          residue_atom_indices( *itr).PushBack( atom_index);
        }
      }

      // create a full atom graph of the molecule
      ConformationGraphConverter::t_AtomGraph atom_graph( ConformationGraphConverter::CreateGraphWithAtoms( *this));

      // for each residue, construct the corresponding molecule and append it to the ensemble
      util::SiPtrVector< const biol::AABase>::const_iterator itr_residues( m_Sequence.Begin());
      util::OwnPtr< graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> > ptr_graph
      (
        &atom_graph,
        false
      );
      for
      (
        storage::Vector< storage::Vector< size_t> >::const_iterator
          itr( residue_atom_indices.Begin()), itr_end( residue_atom_indices.End());
        itr != itr_end;
        ++itr, ++itr_residues
      )
      {
        ensemble.PushBack
        (
          FragmentComplete
          (
            ConformationGraphConverter::CreateAtomsFromGraph
            (
              graph::Subgraph< util::SiPtr< const AtomConformationalInterface>, size_t>( ptr_graph, *itr).ToGraph()
            ),
            ( *itr_residues)->GetIdentification()
          )
        );
      }
      return ensemble;
    }

    //! @brief reconstruct protein model
    //! @param MODEL the original protein model; needed to get the seqres
    //! @note the SSEs present in the model will be ignored entirely
    assemble::ProteinModel AAFragmentComplete::ReconstructProteinModel( const assemble::ProteinModel &MODEL)
    {
      assemble::ProteinModel model;
      storage::Vector< util::ShPtrVector< biol::AABase> > chain_aas( size_t( 127));
      for( auto itr( m_Sequence.Begin()), itr_end( m_Sequence.End()); itr != itr_end; ++itr)
      {
        chain_aas( size_t( ( *itr)->GetChainID())).PushBack( util::CloneToShPtr( **itr));
      }
      for( auto itr( MODEL.GetChains().Begin()), itr_end( MODEL.GetChains().End()); itr != itr_end; ++itr)
      {
        util::ShPtr< assemble::Chain> chain( new assemble::Chain( ( *itr)->GetSequence()));
        const util::ShPtrVector< biol::AABase> &reses( chain_aas( size_t( ( *itr)->GetChainID())));
        // inserts each residue as a new sse with type coil because there are no guarantees that
        // the residues are contiguous, which unfortunately is a requirement for them to be in the AA-sequence that goes
        // into an SSE (granted, we wouldn't need an SSE if we could easily print out an arbitrary set of residues, but
        // that's not currently supported so we have to use a workaround).
        for( auto itr_res( reses.Begin()), itr_res_end( reses.End()); itr_res != itr_res_end; ++itr_res)
        {
          util::ShPtr< assemble::SSE> sse_sp
          (
            new assemble::SSE
            (
              biol::AASequence
              (
                util::ShPtrVector< biol::AABase>( size_t( 1), *itr_res),
                ( *itr)->GetChainID()
              ),
              biol::GetSSTypes().COIL
            )
          );
          chain->Insert( sse_sp);
        }
        model.Insert( chain);
      }
      return model;
    }

    //! @brief get the number of valences
    //! @return the number of valences (=0)
    void AAFragmentComplete::RemoveAtomsUndefinedPos()
    {
      storage::Vector< size_t> defined_ids;
      defined_ids.AllocateMemory( this->GetNumberAtoms());
      auto itr( this->GetAtomsIterator());
      for( size_t i( 0); itr.NotAtEnd(); ++i, ++itr)
      {
        if( itr->GetPosition().IsDefined())
        {
          defined_ids.PushBack( i);
        }
      }
      this->m_Atoms.Reorder( defined_ids);
      m_BiolAtomType.Reorder( defined_ids);
      m_AtomAASequence.Reorder( defined_ids);
      m_AAIndices.Reorder( defined_ids);
      this->ResetCache();
    }

    //! @brief construct this object from a siptr vector of AAs
    //! @param SEQUENCE the sequence of AAs of interest
    //! @param ALLOW_UNDEFINED_POSITIONS true to allow undefined atom positions
    //! Default is to assert and fail if heavy atom position is absent
    void AAFragmentComplete::Construct
    (
      const util::SiPtrVector< const biol::AABase> &AA_SEQUENCE,
      const bool &ALLOW_UNDEFINED_POSITIONS
    )
    {
      m_Sequence.Reset();
      m_Sequence.AllocateMemory( AA_SEQUENCE.GetSize());

      std::string name( "Residues");

      // add all residues with defined type and coordinates to the sequence vector
      // also add up the number of atoms total
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          itr( AA_SEQUENCE.Begin()), itr_end( AA_SEQUENCE.End());
        itr != itr_end;
        ++itr
      )
      {
        if
        (
          itr->IsDefined()
          && ( ALLOW_UNDEFINED_POSITIONS || ( *itr)->GetCenter().IsDefined())
          && ( *itr)->GetType().IsDefined()
        )
        {
          m_Sequence.PushBack( *itr);
          name += " " + ( *itr)->GetType()->GetThreeLetterCode();
        }
      }

      m_AtomAASequence.Reset();

      char current_chain( '*');

      AtomVector< AtomComplete> vector;

      // track whether or not each atom is a backbone atom
      m_AAIndices.Reset();
      m_BiolAtomType.Reset();

      // track string names of all atom types
      storage::Vector< std::string> atom_id_names;

      size_t atoms_so_far( 0), previous_type_c_terminus( 0), aa_index( 0);
      // assemble the AAs
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          itr( m_Sequence.Begin()), itr_end( m_Sequence.End());
        itr != itr_end;
        ++itr, ++aa_index
      )
      {
        const biol::AABase &base( **itr);
        const bool is_n_terminal( current_chain != base.GetChainID());
        current_chain = base.GetChainID();

        // determine whether this AA represents the C-terminus
        util::SiPtrVector< const biol::AABase>::const_iterator itr_next( itr);
        ++itr_next;
        while( itr_next != itr_end && !itr_next->IsDefined())
        {
          ++itr_next;
        }
        const bool is_c_terminal( itr_next == itr_end || ( *itr_next)->GetChainID() != current_chain);

        const FragmentConfigurationShared &fragment( ( *itr)->GetType()->GetFragment( is_c_terminal, is_n_terminal));
        const storage::Vector< size_t>
          &types_map( GetAATypeAtomTypeToIndexMap( is_c_terminal, is_n_terminal)( base.GetType().GetIndex()));
        storage::Vector< sdf::AtomInfo> atom_infos( fragment.GetAtomInfo());
        storage::Vector< char> is_defined_coord( atom_infos.GetSize(), char( 0));
        storage::Vector< std::string> local_names( atom_infos.GetSize());
        storage::Vector< biol::AtomType> local_types( atom_infos.GetSize());
        for
        (
          util::SiPtrVector< const biol::Atom>::const_iterator
            itr_atom( base.GetAtoms().Begin()), itr_atom_end( base.GetAtoms().End());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          // get the coordinates of this atom
          const linal::Vector3D &coordinates( ( *itr_atom)->GetCoordinates());
          // get the index of this atom
          const size_t atom_index( types_map( ( *itr_atom)->GetType().GetIndex()));
          std::string resatom_name
          (
            util::Format().W( 5)( ( *itr)->GetSeqID()) + "_" +
            ( *itr)->GetType()->GetOneLetterCode() + "_" +
            ( *itr)->GetType()->GetThreeLetterCode() + "_" +
            ( *itr_atom)->GetType().GetName()
          );
          if
          (
            !util::IsDefined( atom_index)
            && ( *itr_atom)->GetType()->GetElementType() == GetElementTypes().e_Hydrogen
          )
          {
            BCL_MessageDbg( "Ignoring likely ionized H: " + resatom_name);
            continue;
          }
          BCL_Assert
          (
            util::IsDefined( atom_index),
            "Undefined atom index for type: " + ( *itr_atom)->GetType().GetName()
            + " on residue: " + base.GetType().GetName()
            + " c-term? " + util::Format()( is_c_terminal)
            + " n-term? " + util::Format()( is_n_terminal)
          );
          local_types( atom_index) = ( *itr_atom)->GetType();
          local_names( atom_index) = resatom_name;
          is_defined_coord( atom_index) = '1';
          atom_infos( atom_index).SetCoordinates( coordinates);
        }
        storage::Vector< size_t> defined_atoms;
        defined_atoms.AllocateMemory( atom_infos.GetSize());
        size_t number_undefined_before_c_terminus( 0);
        const size_t c_terminal_id( types_map( biol::GetAtomTypes().C));
        for( size_t i( 0), n_atoms( atom_infos.GetSize()); i < n_atoms; ++i)
        {
          if( is_defined_coord( i))
          {
            defined_atoms.PushBack( i);
          }
          else if( ALLOW_UNDEFINED_POSITIONS && atom_infos( i).GetAtomType()->GetElementType() != GetElementTypes().e_Hydrogen)
          {
            defined_atoms.PushBack( i);
            atom_infos( i).SetCoordinates( linal::Vector3D( util::GetUndefined< double>()));
          }
          else if( !is_c_terminal && i < c_terminal_id)
          {
            ++number_undefined_before_c_terminus;
          }
          BCL_Assert
          (
            is_defined_coord( i)
            || atom_infos( i).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen
            || ALLOW_UNDEFINED_POSITIONS,
            "Missing heavy atom of type " + atom_infos( i).GetAtomType().GetName() + " on residue " + base.GetIdentification()
          );
        }
        const size_t new_c_terminus_id( vector.GetSize() + c_terminal_id - number_undefined_before_c_terminus);
        AtomVector< AtomComplete> atoms( atom_infos, fragment.GetBondInfo());

        if( defined_atoms.GetSize() != atom_infos.GetSize())
        {
          atoms.Reorder( defined_atoms);
          local_names.Reorder( defined_atoms);
          local_types.Reorder( defined_atoms);
        }

        // update the vector that contains a pointer for each atom
        // do not replace storage::Vector< util::SiPtr< const biol::AABase> > with
        // util::SiPtrVector< const biol::AABase>; this would iterate itr_atom_aa, rather than just
        // copying it N times, as we want here
        m_AtomAASequence.Append( storage::Vector< util::SiPtr< const biol::AABase> >( atoms.GetSize(), *itr));
        m_AAIndices.Append( storage::Vector< size_t>( atoms.GetSize(), aa_index));
        atom_id_names.Append( local_names);
        m_BiolAtomType.Append( local_types);

        if( is_n_terminal)
        {
          vector.AddAtomsWithConnectivity( atoms, storage::Vector< sdf::BondInfo>());
        }
        else
        {
          vector.AddAtomsWithConnectivity
          (
            atoms,
            storage::Vector< sdf::BondInfo>
            (
              size_t( 1),
              sdf::BondInfo
              (
                previous_type_c_terminus,
                vector.GetSize(),
                GetConfigurationalBondTypes().e_NonConjugatedSingleBond
              )
            )
          );
        }
        if( !is_c_terminal)
        {
          previous_type_c_terminus = new_c_terminus_id;
        }
        atoms_so_far += fragment.GetNumberAtoms();
      }

      // determine atom and bond types
      AtomsCompleteStandardizer standardizer( vector, name, false);

      // add bond isometry information
      BondIsometryHandler::AddIsometryInformation( vector, true);

      // add stereocenter information
      StereocentersHandler::AddChiralityFromConformation( vector);

      FragmentComplete::operator=( FragmentComplete( vector, name));

      // count the number of valences on each residues, so that we can report the number of missing H on each residue
      // and update m_AtomAASequence accordingly
      util::SiPtrVector< const biol::AABase>::const_iterator itr_atom_aa( m_AtomAASequence.Begin());
      util::SiPtrVector< const biol::AABase> residues_missing_h;
      BCL_Assert
      (
        m_AtomAASequence.GetSize() == GetNumberAtoms(),
        "wrong atom aa sequence size; atom aa sequence size: " + util::Format()( m_AtomAASequence.GetSize())
        + " # atoms: " + util::Format()( GetNumberAtoms())
      );
      storage::Vector< size_t>::const_iterator itr_aa_index( m_AAIndices.Begin());
      storage::Vector< size_t> res_ids_missing_h;
      storage::Vector< biol::AtomType>::const_iterator itr_biol_atom_type( m_BiolAtomType.Begin());
      storage::Vector< biol::AtomType> missing_h_types;
      for
      (
        AtomVector< AtomComplete>::const_iterator
          itr_atom( GetAtomVector().Begin()), itr_atom_end( GetAtomVector().End());
        itr_atom != itr_atom_end;
        ++itr_atom, ++itr_atom_aa, ++itr_aa_index, ++itr_biol_atom_type
      )
      {
        const size_t number_valence_bonds( itr_atom->GetNumberValenceBonds());
        if( number_valence_bonds)
        {
          // do not replace storage::Vector< util::SiPtr< const biol::AABase> > with
          // util::SiPtrVector< const biol::AABase>; this would iterate itr_atom_aa, rather than just
          // copying it N times, as we want here
          residues_missing_h.Append
          (
            storage::Vector< util::SiPtr< const biol::AABase> >( number_valence_bonds, *itr_atom_aa)
          );
          res_ids_missing_h.Append( storage::Vector< size_t>( number_valence_bonds, *itr_aa_index));
          atom_id_names.Append
          (
            storage::Vector< std::string>
            (
              number_valence_bonds,
              std::string
              (
                util::Format().W( 5)( ( *itr_atom_aa)->GetSeqID()) + "_" +
                ( *itr_atom_aa)->GetType()->GetOneLetterCode() + "_" +
                ( *itr_atom_aa)->GetType()->GetThreeLetterCode() + "_" +
                "HXX"
              )
            )
          );

          // determine what type of H is needed
          size_t number_remaining_h( number_valence_bonds);
          for
          (
            storage::Set< biol::AtomType>::const_iterator
              itr_type( ( *itr_biol_atom_type)->GetConnections().Begin()),
              itr_type_end( ( *itr_biol_atom_type)->GetConnections().End());
            itr_type != itr_type_end;
            ++itr_type
          )
          {
            if( ( *itr_type)->GetElementType() == GetElementTypes().e_Hydrogen)
            {
              missing_h_types.PushBack( *itr_type);
              if( !--number_remaining_h)
              {
                break;
              }
            }
          }
          // if none is found (which is the case for various Hs of MTSL for which we lack the atom type), then just add HH
          if( number_remaining_h)
          {
            missing_h_types.Append( storage::Vector< biol::AtomType>( number_remaining_h, biol::GetAtomTypes().HH));
          }
          BCL_MessageDbg
          (
            "Residue " + ( *itr_atom_aa)->GetIdentification()
            + " missing " + util::Format()( number_valence_bonds) + " H on a "
            + itr_atom->GetAtomType().GetName()
          );
        }
      }
      m_AtomAASequence.Append( residues_missing_h);
      m_AAIndices.Append( res_ids_missing_h);
      m_BiolAtomType.Append( missing_h_types);

      std::string joined_names;
      for( size_t i( 0), n_atoms( atom_id_names.GetSize()); i < n_atoms; ++i)
      {
        joined_names += atom_id_names( i) + ' ';
      }
      StoreProperty( "SeqResAtomName", joined_names);

      this->SaturateWithH();
    }

  } // namespace chemistry
} // namespace bcl
