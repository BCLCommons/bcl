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
#include "descriptor/bcl_descriptor_molecule_halogenated_aromatic_rings.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeHalogenatedAromaticRings::MoleculeHalogenatedAromaticRings
    (
    ) :
      m_CountMethod( MoleculeHalogenatedAromaticRings::AromaticRingHalogens::e_Total)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief constructor
    //! @brief COUNT_METHOD how to count halogens
    MoleculeHalogenatedAromaticRings::MoleculeHalogenatedAromaticRings
    (
      const MoleculeHalogenatedAromaticRings::AromaticRingHalogens &COUNT_METHOD
    ) :
      m_CountMethod( COUNT_METHOD)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief constructor with halogen types
    //! @brief COUNT_METHOD how to count halogens
    //! @brief ALLOWED_HALOGENS the halogens included in counting
    MoleculeHalogenatedAromaticRings::MoleculeHalogenatedAromaticRings
    (
      const MoleculeHalogenatedAromaticRings::AromaticRingHalogens &COUNT_METHOD,
      const storage::Vector< chemistry::AtomType> &ALLOWED_HALOGENS
    ) :
      m_CountMethod( COUNT_METHOD),
      m_AllowedHalogens( ALLOWED_HALOGENS)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeHalogenatedAromaticRings
    MoleculeHalogenatedAromaticRings *MoleculeHalogenatedAromaticRings::Clone() const
    {
      return new MoleculeHalogenatedAromaticRings( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeHalogenatedAromaticRings::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeHalogenatedAromaticRings::GetAlias() const
    {
      static const std::string s_total( "NAromaticRingHalogensTotal");
      static const std::string s_max( "NAromaticRingHalogensMaxFragment");
      return m_CountMethod == e_Max ? s_max : s_total;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeHalogenatedAromaticRings::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // Get our molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      chemistry::FragmentComplete molecule( *current_mol);

      // count the total number of halogen atoms on aromatic rings
      if( m_CountMethod == e_Total)
      {
        STORAGE = float( CountAromaticHalogens( molecule));
      }
      else if( m_CountMethod == e_Max)
      {
        STORAGE = float( CountMaxAromaticHalogensSingleFragment( molecule));
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the max number of halogen atoms on a ring in a molecule
    //! @param MOLECULE the molecule of interest
    size_t MoleculeHalogenatedAromaticRings::CountMaxAromaticHalogensSingleFragment( const chemistry::FragmentComplete &MOLECULE) const
    {
      // generate a molecule splitter
      chemistry::FragmentSplitRigid splitter;
      chemistry::FragmentEnsemble rigid_fragments( splitter( MOLECULE));

      // see if fragment is a ring
      size_t max_n_halogens( 0);
      for
      (
          auto frag_itr( rigid_fragments.Begin()), frag_itr_end( rigid_fragments.End());
          frag_itr != frag_itr_end;
          ++frag_itr
      )
      {
        // start at 0 for each fragment
        size_t n_halogens( 0);

        // irrelevant if not an aromatic ring
        if( GetCheminfoProperties().calc_NAromaticRings->SumOverObject( *frag_itr)( 0))
        {
          // default to all useful halogens
          if( m_AllowedHalogens.GetSize() == size_t( 4))
          {
            n_halogens += GetCheminfoProperties().calc_IsHalogen->SumOverObject( *frag_itr)( 0);
          }
          // only tally select halogens
          else
          {
            if( m_AllowedHalogens.Find( chemistry::GetAtomTypes().F_SP2P2P2) < m_AllowedHalogens.GetSize())
            {
              n_halogens += GetCheminfoProperties().calc_IsF->SumOverObject( *frag_itr)( 0);
            }
            if( m_AllowedHalogens.Find( chemistry::GetAtomTypes().Cl_SP2P2P2) < m_AllowedHalogens.GetSize())
            {
              n_halogens += GetCheminfoProperties().calc_IsCl->SumOverObject( *frag_itr)( 0);
            }
            if( m_AllowedHalogens.Find( chemistry::GetAtomTypes().Br_SP2P2P2) < m_AllowedHalogens.GetSize())
            {
              n_halogens += GetCheminfoProperties().calc_IsBr->SumOverObject( *frag_itr)( 0);
            }
            if( m_AllowedHalogens.Find( chemistry::GetAtomTypes().I_SP2P2P2) < m_AllowedHalogens.GetSize())
            {
              n_halogens += GetCheminfoProperties().calc_IsI->SumOverObject( *frag_itr)( 0);
            }
          }
          // set new max value
          if( n_halogens > max_n_halogens)
          {
            max_n_halogens = n_halogens;
          }
        }
      }
      return max_n_halogens;
    }

    //! @brief calculate the total number of non-aromatic ring Cl, Br, I halogen atoms
    //! @param MOLECULE the molecule of interest
    size_t MoleculeHalogenatedAromaticRings::CountAromaticHalogens( const chemistry::FragmentComplete &MOLECULE) const
    {
      // generate a molecule splitter
      chemistry::FragmentSplitRigid splitter;
      chemistry::FragmentEnsemble rigid_fragments( splitter( MOLECULE));

      // see if fragment is an aromatic ring
      size_t n_halogens( 0);
      for
      (
          auto frag_itr( rigid_fragments.Begin()), frag_itr_end( rigid_fragments.End());
          frag_itr != frag_itr_end;
          ++frag_itr
      )
      {
        // find aromatic rings
        if( GetCheminfoProperties().calc_NAromaticRings->SumOverObject( *frag_itr)( 0))
        {
          // default to all halogens
          if( m_AllowedHalogens.GetSize() == size_t( 4))
          {
            n_halogens += GetCheminfoProperties().calc_IsHalogen->SumOverObject( *frag_itr)( 0);
          }
          else
          {
            if( m_AllowedHalogens.Find( chemistry::GetAtomTypes().F_SP2P2P2) < m_AllowedHalogens.GetSize())
            {
              n_halogens += GetCheminfoProperties().calc_IsF->SumOverObject( *frag_itr)( 0);
            }
            if( m_AllowedHalogens.Find( chemistry::GetAtomTypes().Cl_SP2P2P2) < m_AllowedHalogens.GetSize())
            {
              n_halogens += GetCheminfoProperties().calc_IsCl->SumOverObject( *frag_itr)( 0);
            }
            if( m_AllowedHalogens.Find( chemistry::GetAtomTypes().Br_SP2P2P2) < m_AllowedHalogens.GetSize())
            {
              n_halogens += GetCheminfoProperties().calc_IsBr->SumOverObject( *frag_itr)( 0);
            }
            if( m_AllowedHalogens.Find( chemistry::GetAtomTypes().I_SP2P2P2) < m_AllowedHalogens.GetSize())
            {
              n_halogens += GetCheminfoProperties().calc_IsI->SumOverObject( *frag_itr)( 0);
            }
          }
        }
      }
      return n_halogens;
    }

    //! @brief set the halogens to be counted
    void MoleculeHalogenatedAromaticRings::SetAllowedHalogens( const storage::Vector< chemistry::AtomType> &ALLOWED_HALOGENS)
    {
      m_AllowedHalogens = ALLOWED_HALOGENS;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool MoleculeHalogenatedAromaticRings::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // static initialization check
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }

      // read in allowed halogens
      if( m_AllowedHalogensString.size())
      {
        // parse input
        const storage::Vector< std::string> allowed_halogens
        (
          util::SplitString( util::TrimString( m_AllowedHalogensString), " \t\n\r,")
        );

        // stupid check to add only the correct halogens
        for( size_t h_i( 0), h_sz( allowed_halogens.GetSize()); h_i < h_sz; ++h_i)
        {
          // Fluorine
          if( allowed_halogens( h_i) == "F")
          {
            m_AllowedHalogens.PushBack( chemistry::GetAtomTypes().F_SP2P2P2);
          }
          // Chlorine
          if( allowed_halogens( h_i) == "Cl")
          {
            m_AllowedHalogens.PushBack( chemistry::GetAtomTypes().Cl_SP2P2P2);
          }
          // Bromine
          if( allowed_halogens( h_i) == "Br")
          {
            m_AllowedHalogens.PushBack( chemistry::GetAtomTypes().Br_SP2P2P2);
          }
          // Iodine
          if( allowed_halogens( h_i) == "I")
          {
            m_AllowedHalogens.PushBack( chemistry::GetAtomTypes().I_SP2P2P2);
          }
        }
      }

      // done
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeHalogenatedAromaticRings::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Determines the sum of the number of halogen substituents on aromatic rings in a molecule. "
        "Alternatively, identifies the aromatic ring with the most number of halogens and returns that count."
      );
      parameters.AddInitializer
      (
        "allowed_halogens",
        "halogen types that will be counted",
        io::Serialization::GetAgent( &m_AllowedHalogensString),
        "F Cl Br I"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
