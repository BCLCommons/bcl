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
#include "chemistry/bcl_chemistry_mutate_bond_lengths.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_clash_score.h"
#include "chemistry/bcl_chemistry_bond_angle_assignment.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_mutate_bond_angles.h"
#include "chemistry/bcl_chemistry_mutate_chirality.h"
#include "chemistry/bcl_chemistry_mutate_clash_resolver.h"
#include "chemistry/bcl_chemistry_mutate_dihedrals_interface.h"
#include "chemistry/bcl_chemistry_possible_atom_types_for_atom.h"
#include "chemistry/bcl_chemistry_ring_fragment_map.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_mutate_result.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

#define BCL_PROFILE_MutateBondLengths
#ifdef BCL_PROFILE_MutateBondLengths
#include "util/bcl_util_stopwatch.h"
#endif

namespace bcl
{
  namespace chemistry
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateBondLengths::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateBondLengths())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new MutateBondLengths
    MutateBondLengths *MutateBondLengths::Clone() const
    {
      return new MutateBondLengths( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateBondLengths::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateBondLengths::GetScheme() const
    {
      static const std::string s_name( "MinimizeBondLengths");
      return s_name;
    }

    //! @brief get access to the stored molecule info
    //! @return the MolBondInfo object
    const MutateBondLengths::MolBondInfo &MutateBondLengths::GetInfo() const
    {
      return m_Info;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking a conformation and returns a mutated conformation
    //! @param MOLECULE conformation of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< FragmentComplete> MutateBondLengths::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      // get info on starting molecule
      FragmentComplete molecule( MOLECULE);
      if( !m_Initialized)
      {
        Initialize( molecule);
      }

      // bail if our bond lengths are all near ideal values
      if( !util::IsDefined( m_Info.m_WorstBondIndex))
      {
        // return empty; mutate is skipped
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // start the timer
      #ifdef BCL_PROFILE_MutateBondLengths
      static util::Stopwatch s_timer( "Minimize bond lengths", util::Message::e_Standard, true);
      s_timer.Start();
      #endif

      // mobile atoms
      if( !m_MobileAtoms.GetSize())
      {
        m_MobileAtoms = GetMobileAtoms( molecule);
      }
      // resolve any nan value sin atom types
      for( size_t i( 0); i < m_MaxCycles; ++i)
      {
        UpdateAtomTypes( m_Info, m_MobileAtoms);
        UpdateInfo( m_Info, m_MobileAtoms);
        molecule = FragmentComplete( m_Info.m_Atoms, molecule.GetName());
        molecule.StandardizeBondLengths( m_MobileAtoms);
      }

      // stop timer
      #ifdef BCL_PROFILE_MutateBondLengths
      s_timer.Stop();
      #endif

      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>( new FragmentComplete( molecule.GetAtomVector(), MOLECULE.GetName())), *this);

//      BCL_MessageDbg( "Done...Clash: " + util::Format()( std::min( clash_score, old_cs)));
//      if( old_cs < clash_score)
//      {
//        return math::MutateResult< FragmentComplete>( best_mol_sp, *this);
//      }
//
//      return math::MutateResult< FragmentComplete>( current_sp, *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize info object with molecule
    void MutateBondLengths::Initialize( const FragmentComplete &MOLECULE) const
    {
      m_Info = MutateBondLengths::MolBondInfo( MOLECULE);
      m_Initialized = true;
    }

    //! @brief setup this mutate to handle a new molecule
    storage::Vector< size_t> MutateBondLengths::GetMobileAtoms
    (
      const FragmentComplete &CONF
    ) const
    {
      // Get indices that we are allowed to sample
      linal::Vector< size_t> sample_by_parts( CONF.GetStoredProperties().GetMDLPropertyAsVector( "SampleByParts"));
      storage::Set< size_t> parts_to_sample;

      // get the unique atom indices from the MDL property
      if( sample_by_parts.GetSize())
      {
        for
        (
            linal::Vector< size_t>::const_iterator itr_indices( sample_by_parts.Begin()), itr_indices_end( sample_by_parts.End());
            itr_indices != itr_indices_end;
            ++itr_indices
        )
        {
          parts_to_sample.Insert( *itr_indices);
        }
      }
      else
      {
        for( size_t i( 0); i < CONF.GetSize(); ++i)
        {
          parts_to_sample.Insert( i);
        }
      }

      // get all atom indices
      return storage::Vector< size_t>( parts_to_sample.Begin(), parts_to_sample.End());
    }

    //! @brief perturb the molecule to update bond lengths
    //! @param INFO bond length accounting info of current molecule
    //! @param MOBILE_ATOMS atoms that can be moved during perturbation
    void MutateBondLengths::UpdateBondLengths
    (
      MutateBondLengths::MolBondInfo &INFO,
      const storage::Vector< size_t> MOBILE_ATOMS
    ) const
    {
      // nudge an atom from each pair along its bond vector
      size_t bond_index( 0);
      for
      (
          auto bond_itr( INFO.m_BondInfo.Begin()), bond_itr_end( INFO.m_BondInfo.End());
          bond_itr != bond_itr_end;
          ++bond_itr, ++bond_index
      )
      {
        // reference the atoms in the bond
        const AtomComplete &atom_a( INFO.m_Atoms( bond_itr->GetAtomIndexLow()));
        const AtomComplete &atom_b( INFO.m_Atoms( bond_itr->GetAtomIndexHigh()));

        // ignore hydrogen atom bonds
        if
        (
            atom_a.GetElementType() == GetElementTypes().e_Hydrogen ||
            atom_b.GetElementType() == GetElementTypes().e_Hydrogen ||
            !INFO.m_BondLengthsDiff( bond_index)
        )
        {
          continue;
        }

        BCL_MessageStd
        (
          "Bond between atoms " +
          util::Format()( INFO.m_Atoms.GetAtomIndex( atom_a)) +
          " and " +
          util::Format()( INFO.m_Atoms.GetAtomIndex( atom_b))
        );
        BCL_Debug( INFO.m_BondLengthsCurrent( bond_index));
        BCL_Debug( INFO.m_BondLengthsEq( bond_index));
        BCL_Debug( INFO.m_BondLengthsDiff( bond_index));

        // create bond vector
        linal::Vector3D old_bond_vec( atom_a.GetPosition() - atom_b.GetPosition());
        linal::Vector3D bond_vec( old_bond_vec);
        BCL_Debug( bond_vec);
        double norm( bond_vec.Norm());
        BCL_Debug( norm);
        BCL_Debug( m_StepSize / norm - 1.0);
        bond_vec *= m_StepSize / norm - 1.0;
        BCL_Debug( bond_vec);

        // adjust the position of a mobile atom in the current pair
        size_t atom_index( util::GetUndefinedSize_t());
        if( MOBILE_ATOMS.Find( bond_itr->GetAtomIndexLow()) < MOBILE_ATOMS.GetSize())
        {
          atom_index = bond_itr->GetAtomIndexLow();
        }
        else if( MOBILE_ATOMS.Find( bond_itr->GetAtomIndexHigh()) < MOBILE_ATOMS.GetSize())
        {
          atom_index = bond_itr->GetAtomIndexHigh();
        }

        if( util::IsDefined( atom_index))
        {
          BCL_MessageStd( "Position before update: " + util::Format()( INFO.m_Atoms( atom_index).GetPosition()));
          linal::Vector3D pos( INFO.m_Atoms( atom_index).GetPosition());
          linal::Vector3D new_pos( pos + ( bond_vec + old_bond_vec) * -1.0);
          INFO.m_Atoms( atom_index).SetPosition( new_pos);
          BCL_MessageStd( "Position after update: " + util::Format()( INFO.m_Atoms( atom_index).GetPosition()));
        }

        // adjust the positions of all mobile atoms
        for
        (
            auto atom_itr( MOBILE_ATOMS.Begin()), atom_itr_end( MOBILE_ATOMS.End());
            atom_itr != atom_itr_end;
            ++atom_itr
        )
        {
          linal::Vector3D pos( INFO.m_Atoms( *atom_itr).GetPosition());
          linal::Vector3D new_pos( pos + ( bond_vec + old_bond_vec) * -1.0);
          INFO.m_Atoms( *atom_itr).SetPosition( new_pos);
        }

      }

      FragmentComplete update_frag( INFO.m_Atoms, "");

      io::OFStream debug_out;
      io::File::MustOpenOFStream( debug_out, "REACTION.updates.sdf", std::ios::app);
      update_frag.WriteMDL( debug_out);
      io::File::CloseClearFStream( debug_out);
    }

    //! @brief update atom types of select atoms
    //! @param INFO bond types accounting info of current molecule
    //! @param MOBILE_ATOMS atoms whose bond type can be altered
    void MutateBondLengths::UpdateAtomTypes
    (
      MolBondInfo &INFO,
      const storage::Vector< size_t> MOBILE_ATOMS
    ) const
    {
      // fix undefined bond lengths
      for( size_t i( 0); i < INFO.m_NumberBonds; ++i)
      {
        // get atom reference
        auto &atom_a( INFO.m_Atoms( INFO.m_BondInfo( i).GetAtomIndexLow()));
        auto &atom_b( INFO.m_Atoms( INFO.m_BondInfo( i).GetAtomIndexHigh()));

        // equilibrium bond length
        double bond_length
        (
          BondLengths::GetBondLength
          (
            atom_a.GetAtomType(),
            INFO.m_BondInfo( i).GetConfigurationalBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic),
            atom_b.GetAtomType()
          )
        );

        // if the bond length is undefined then find the offending atom
        if( !util::IsDefined( bond_length))
        {
          double bond_a
          (
            BondLengths::GetCovalentRadius
            (
              atom_a.GetAtomType(),
              INFO.m_BondInfo( i).GetConfigurationalBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic)
            )
          );
          if( !util::IsDefined( bond_a))
          {
            AtomType new_atom_type( FindStableAtomType( atom_a));
            atom_a.SetAtomType( new_atom_type);
          }
          double bond_b
          (
            BondLengths::GetCovalentRadius
            (
              atom_b.GetAtomType(),
              INFO.m_BondInfo( i).GetConfigurationalBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic)
            )
          );
          if( !util::IsDefined( bond_b))
          {
            AtomType new_atom_type( FindStableAtomType( atom_b));
            atom_b.SetAtomType( new_atom_type);
          }
        }
      }
    }

    //! @brief get a new stable atom type
    AtomType MutateBondLengths::FindStableAtomType
    (
      const util::SiPtr< const AtomConformationalInterface> &ATOM
    ) const
    {
      // starting atom type
      AtomType atom_type( ATOM->GetAtomType());

      // aromatic?
      bool picked_atom_aromatic( ATOM->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, 1));

      // find number of bonds being made by starting atom
      size_t n_h( ATOM->GetNumberCovalentlyBoundHydrogens());
      size_t n_nonh_bonds( ATOM->GetBonds().GetSize() - n_h);
      for( size_t n_current_bonds( n_nonh_bonds); n_current_bonds <= 6; ++n_current_bonds)
      {
        size_t n_nonh_e( ATOM->GetAtomType()->GetNumberElectronsInBonds() - n_h + ( n_current_bonds - n_nonh_bonds));
        PossibleAtomTypesForAtom available_atom_types
        (
          ATOM->GetElementType(),
          n_nonh_e,
          n_current_bonds,
          ATOM->GetCharge(),
          picked_atom_aromatic
        );
        // prefer original charge
        if( available_atom_types.GetNumberPossibleTypes() && available_atom_types.GetMostStableType()->IsGasteigerAtomType())
        {
          if( available_atom_types.GetAlternateTypeWithCharge( ATOM->GetCharge()).IsDefined())
          {
            atom_type = available_atom_types.GetAlternateTypeWithCharge( ATOM->GetCharge());
            break;
          }
          atom_type = available_atom_types.GetMostStableType();
          break;
        }
      }
      return atom_type;
    }

    //! @brief update the MolBondInfo object
    //! @param INFO bond length accounting info of current molecule
    //! @param MOBILE_ATOMS atoms that can be moved during perturbation
    void MutateBondLengths::UpdateInfo
    (
      MutateBondLengths::MolBondInfo &INFO,
      const storage::Vector< size_t> MOBILE_ATOMS
    ) const
    {
      // set values for each bond
      for( size_t i( 0); i < INFO.m_NumberBonds; ++i)
      {
        // get atom reference
        auto &atom_a( INFO.m_Atoms( INFO.m_BondInfo( i).GetAtomIndexLow()));
        auto &atom_b( INFO.m_Atoms( INFO.m_BondInfo( i).GetAtomIndexHigh()));

        // equilibrium bond length
        INFO.m_BondLengthsEq( i) = BondLengths::GetBondLength
        (
          atom_a.GetAtomType(),
          INFO.m_BondInfo( i).GetConfigurationalBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic),
          atom_b.GetAtomType()
        );

        // current bond length
        linal::Vector3D bond_vec( atom_a.GetPosition() - atom_b.GetPosition());
        INFO.m_BondLengthsCurrent( i) = bond_vec.Norm();

        // difference from equilibrium accounting for tolerance
        INFO.m_BondLengthsDiff( i) =
            math::Absolute( INFO.m_BondLengthsCurrent( i) - INFO.m_BondLengthsEq( i)) > INFO.m_BondLengthsTolerance( i) ?
                INFO.m_BondLengthsCurrent( i) - INFO.m_BondLengthsEq( i) :
                0.0;
      }

      // get the worst bond index
      double greatest_deviation( 0.0);
      for( size_t i( 0); i < INFO.m_NumberBonds; ++i)
      {
        if( math::Absolute( INFO.m_BondLengthsDiff( i)) > greatest_deviation)
        {
          greatest_deviation = math::Absolute( INFO.m_BondLengthsDiff( i));
          INFO.m_WorstBondIndex = i;
        }
      }

      // if all of the bonds are good then set worst undefined
      if( !greatest_deviation)
      {
        INFO.m_WorstBondIndex = util::GetUndefinedSize_t();
      }
    }

    //! @brief check if all bond lengths meet convergence criteria
    //! @param INFO bond length accounting info of the current molecule
    //! @param MOBILE_ATOMS atoms contributing bonds that can be perturbed
    //! @return true if all bond lengths are within the convergence limit
    bool MutateBondLengths::CheckConvergence
    (
      const MutateBondLengths::MolBondInfo &INFO,
      const storage::Vector< size_t> MOBILE_ATOMS
    ) const
    {
      for( size_t i( 0); i < INFO.m_NumberBonds; ++i)
      {
        // skip if neither bonded atom is a mobile atom
        if
        (
            MOBILE_ATOMS.Find( INFO.m_BondInfo( i).GetAtomIndexLow()) == MOBILE_ATOMS.GetSize() &&
            MOBILE_ATOMS.Find( INFO.m_BondInfo( i).GetAtomIndexHigh()) == MOBILE_ATOMS.GetSize()
        )
        {
          continue;
        }
        // check convergence criteria
        if( math::Absolute( INFO.m_BondLengthsDiff( i)) > m_ConvergenceLimit)
        {
          return false;
        }
      }
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateBondLengths::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription( "Resolves clashes");
      serial.AddInitializer
      (
        "max cycles",
        "maximum number of times clashes are redected and removed",
        io::Serialization::GetAgent( &m_MaxCycles),
        "5"
      );
      return serial;
    }

  } // namespace chemistry
} // namespace bcl
