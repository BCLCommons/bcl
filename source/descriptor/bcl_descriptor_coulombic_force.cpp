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
#include "descriptor/bcl_descriptor_coulombic_force.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_identity.h"
#include "math/bcl_math_trigonometric_transition.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> CoulombicForce::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new CoulombicForce( false)
      )
    );

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> CoulombicForce::s_DihedralInstance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new CoulombicForce( true)
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from whether to focus only on dihedral interactions
    //! @param DIHEDRAL_INTERACTIONS_ONLY if true, focus only on dihedral interactions
    CoulombicForce::CoulombicForce( bool DIHEDRAL_INTERACTIONS_ONLY) :
      m_ExcludeNeighbors( DIHEDRAL_INTERACTIONS_ONLY),
      m_ExcludeNeighborsOfNeighbors( DIHEDRAL_INTERACTIONS_ONLY),
      m_DistanceCutoff( util::GetUndefined< double>())
    {
    }

    //! @brief constructor from number of steps, and mapped atom property
    CoulombicForce::CoulombicForce
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const double &DISTANCE_CUTOFF
    ) :
      m_Charge( ATOM_PROPERTY),
      m_ExcludeNeighbors( false),
      m_ExcludeNeighborsOfNeighbors( false),
      m_DistanceCutoff( DISTANCE_CUTOFF)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new CoulombicForce
    CoulombicForce *CoulombicForce::Clone() const
    {
      return new CoulombicForce( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &CoulombicForce::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &CoulombicForce::GetAlias() const
    {
      static const std::string s_name( "CoulombicForce"), s_dihedral_name( "RotamerCoulombicForce");
      return m_ExcludeNeighborsOfNeighbors ? s_dihedral_name : s_name;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &CoulombicForce::GetAtomProperty() const
    {
      return m_Charge;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > CoulombicForce::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_Charge,
          &m_Charge + 1
        );
    }

    //! @brief get intermolecular force
    //! @param MOL_A, MOL_B molecules of interest
    double CoulombicForce::GetIntermolecularForce
    (
      const chemistry::ConformationInterface &MOL_A,
      const chemistry::ConformationInterface &MOL_B
    )
    {
      const linal::Vector< float> charges_a( m_Charge->CollectValuesOnEachElementOfObject( MOL_A));
      const linal::Vector< float> charges_b( m_Charge->CollectValuesOnEachElementOfObject( MOL_B));
      float force_sum( 0.0);
      auto itr_charge_a( charges_a.Begin());
      const double dist_cutoff_sqr
      (
        util::IsDefined( m_DistanceCutoff) ? math::Sqr( m_DistanceCutoff) : math::GetHighestBoundedValue< double>()
      );
      math::TrigonometricTransition trig;
      math::Identity< double> identity_func;
      if( util::IsDefined( m_DistanceCutoff))
      {
        trig = math::TrigonometricTransition( math::Sqr( m_DistanceCutoff - 1.0), math::Sqr( m_DistanceCutoff), 1.0, 0.0);
      }
      const math::FunctionInterfaceSerializable< double, double> &transition
      (
        util::IsDefined( m_DistanceCutoff)
        ? ( const math::FunctionInterfaceSerializable< double, double> &)trig
        : ( const math::FunctionInterfaceSerializable< double, double> &)identity_func
      );
      for
      (
        auto itr_atoms_a( MOL_A.GetIterator());
        itr_atoms_a.NotAtEnd();
        ++itr_atoms_a, ++itr_charge_a
      )
      {
        if( !itr_atoms_a->GetPosition().IsDefined())
        {
          continue;
        }
        auto itr_charge_b( charges_b.Begin());
        for
        (
          auto itr_atoms_b( MOL_B.GetIterator());
          itr_atoms_b.NotAtEnd();
          ++itr_atoms_b, ++itr_charge_b
        )
        {
          if( !itr_atoms_b->GetPosition().IsDefined())
          {
            continue;
          }

          // store distance between both atoms
          const float sqr_distance( linal::SquareDistance( itr_atoms_a->GetPosition(), itr_atoms_b->GetPosition()));
          if( sqr_distance < 0.1 || sqr_distance > dist_cutoff_sqr)
          {
            // overlapping atom
            continue;
          }
          force_sum += transition( sqr_distance) * *itr_charge_b * *itr_charge_a / sqr_distance;
        }
      }
      return force_sum;
    }

    //! @brief get absolute sum of intermolecular forces (to weight molecular interaction signficance)
    //! @param MOL_A, MOL_B molecules of interest
    double CoulombicForce::GetSumAbsIntermolecularForce
    (
      const chemistry::ConformationInterface &MOL_A,
      const chemistry::ConformationInterface &MOL_B
    )
    {
      const linal::Vector< float> charges_a( m_Charge->CollectValuesOnEachElementOfObject( MOL_A));
      const linal::Vector< float> charges_b( m_Charge->CollectValuesOnEachElementOfObject( MOL_B));
      float force_sum( 0.0);
      auto itr_charge_a( charges_a.Begin());
      const double dist_cutoff_sqr
      (
        util::IsDefined( m_DistanceCutoff) ? math::Sqr( m_DistanceCutoff) : math::GetHighestBoundedValue< double>()
      );
      math::TrigonometricTransition trig;
      math::Identity< double> identity_func;
      if( util::IsDefined( m_DistanceCutoff))
      {
        trig = math::TrigonometricTransition( math::Sqr( m_DistanceCutoff - 1.0), math::Sqr( m_DistanceCutoff), 1.0, 0.0);
      }
      const math::FunctionInterfaceSerializable< double, double> &transition
      (
        util::IsDefined( m_DistanceCutoff)
        ? ( const math::FunctionInterfaceSerializable< double, double> &)trig
        : ( const math::FunctionInterfaceSerializable< double, double> &)identity_func
      );
      for
      (
        auto itr_atoms_a( MOL_A.GetIterator());
        itr_atoms_a.NotAtEnd();
        ++itr_atoms_a, ++itr_charge_a
      )
      {
        if( !itr_atoms_a->GetPosition().IsDefined())
        {
          continue;
        }
        auto itr_charge_b( charges_b.Begin());
        for
        (
          auto itr_atoms_b( MOL_B.GetIterator());
          itr_atoms_b.NotAtEnd();
          ++itr_atoms_b, ++itr_charge_b
        )
        {
          if( !itr_atoms_b->GetPosition().IsDefined())
          {
            continue;
          }

          // store distance between both atoms
          const float sqr_distance( linal::SquareDistance( itr_atoms_a->GetPosition(), itr_atoms_b->GetPosition()));
          if( sqr_distance < 0.1 || sqr_distance > dist_cutoff_sqr)
          {
            // overlapping atom
            continue;
          }
          force_sum += transition( sqr_distance) * math::Absolute( *itr_charge_b * *itr_charge_a) / sqr_distance;
        }
      }
      return force_sum;
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void CoulombicForce::Calculate( linal::VectorReference< float> &STORAGE)
    {
      const linal::Vector< float> charges( m_Charge->CollectValuesOnEachElementOfObject( *( this->GetCurrentObject())));

      float force_sum( 0.0);
      // iterate over all possible pairs of atoms
      // iterate properties and surface areas simultaneously
      linal::Vector< float>::const_iterator itr_charge_a( charges.Begin());
      util::SiPtr< const chemistry::ConformationInterface> conformation( this->GetCurrentObject());

      const double dist_cutoff_sqr
      (
        util::IsDefined( m_DistanceCutoff) ? math::Sqr( m_DistanceCutoff) : math::GetHighestBoundedValue< double>()
      );
      math::TrigonometricTransition trig;
      math::Identity< double> identity_func;
      if( util::IsDefined( m_DistanceCutoff))
      {
        trig = math::TrigonometricTransition( math::Sqr( m_DistanceCutoff - 1.0), math::Sqr( m_DistanceCutoff), 1.0, 0.0);
      }
      const math::FunctionInterfaceSerializable< double, double> &transition
      (
        util::IsDefined( m_DistanceCutoff)
        ? ( const math::FunctionInterfaceSerializable< double, double> &)trig
        : ( const math::FunctionInterfaceSerializable< double, double> &)identity_func
      );

      storage::List< storage::Vector< size_t> > ring_systems;
      util::SiPtrVector< const storage::Vector< size_t> > atom_to_ring_system;
      if( m_ExcludeNeighborsOfNeighbors)
      {
        chemistry::FragmentSplitRings ring_splitter( true, size_t( 3));
        chemistry::ConformationGraphConverter::t_AtomGraph atom_graph
        (
          chemistry::ConformationGraphConverter::CreateGraphWithAtoms( *conformation)
        );
        ring_systems = ring_splitter.GetComponentVertices( *conformation, atom_graph);
        atom_to_ring_system.Resize( conformation->GetSize());
        for
        (
          storage::List< storage::Vector< size_t> >::const_iterator
            itr( ring_systems.Begin()), itr_end( ring_systems.End());
          itr != itr_end;
          ++itr
        )
        {
          util::SiPtr< const storage::Vector< size_t> > siptr_ring( *itr);
          for
          (
            storage::Vector< size_t>::const_iterator itr_atom_id( itr->Begin()), itr_atom_id_end( itr->End());
            itr_atom_id != itr_atom_id_end;
            ++itr_atom_id
          )
          {
            atom_to_ring_system( *itr_atom_id) = siptr_ring;
          }
        }
      }
      std::string is_excluded( conformation->GetSize(), '0');
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface>
          itr_atoms_a( conformation->GetIterator());
        itr_atoms_a.NotAtEnd();
        ++itr_atoms_a, ++itr_charge_a
      )
      {
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_b( itr_atoms_a);
        linal::Vector< float>::const_iterator itr_charge_b( itr_charge_a + 1);

        // populate is_neighbor with 1s for each atom that should be excluded
        // note that for performance, we avoid assigning 0 to all positions, since this is done in the force calculation
        // component below for the positions to the right of itr_atoms_a. The positions to the left of itr_atoms_a in
        // is_excluded are not used
        if( m_ExcludeNeighbors)
        {
          for
          (
            storage::Vector< chemistry::BondConformational>::const_iterator
              itr_bonds_a( itr_atoms_a->GetBonds().Begin()), itr_bonds_a_end( itr_atoms_a->GetBonds().End());
            itr_bonds_a != itr_bonds_a_end;
            ++itr_bonds_a
          )
          {
            is_excluded[ conformation->GetAtomIndex( itr_bonds_a->GetTargetAtom())] = '1';
            if( m_ExcludeNeighborsOfNeighbors)
            {
              for
              (
                storage::Vector< chemistry::BondConformational>::const_iterator
                  itr_bonds_b( itr_bonds_a->GetTargetAtom().GetBonds().Begin()),
                  itr_bonds_b_end( itr_bonds_a->GetTargetAtom().GetBonds().End());
                itr_bonds_b != itr_bonds_b_end;
                ++itr_bonds_b
              )
              {
                is_excluded[ conformation->GetAtomIndex( itr_bonds_b->GetTargetAtom())] = '1';
              }
            }
          }
          const size_t atom_a_id( itr_atoms_a.GetPosition());
          if( m_ExcludeNeighborsOfNeighbors)
          {
            if( itr_atoms_a->GetBonds().GetSize() == size_t( 1))
            {
              atom_to_ring_system( atom_a_id)
                = atom_to_ring_system( conformation->GetAtomIndex( itr_atoms_a->GetBonds().Begin()->GetTargetAtom()));
            }
            if( atom_to_ring_system( atom_a_id).IsDefined())
            {
              for
              (
                storage::Vector< size_t>::const_iterator
                  itr_atom_id( atom_to_ring_system( atom_a_id)->Begin()),
                  itr_atom_id_end( atom_to_ring_system( atom_a_id)->End());
                itr_atom_id != itr_atom_id_end;
                ++itr_atom_id
              )
              {
                is_excluded[ *itr_atom_id] = '1';
              }
            }
          }
        }

        std::string::iterator itr_excluded( is_excluded.begin() + itr_atoms_a.GetPosition() + 1);
        for( ++itr_atoms_b; itr_atoms_b.NotAtEnd(); ++itr_atoms_b, ++itr_charge_b, ++itr_excluded)
        {
          if( *itr_excluded == '1')
          {
            *itr_excluded = '0';
            continue;
          }
          // store distance between both atoms
          const float sqr_distance( linal::SquareDistance( itr_atoms_a->GetPosition(), itr_atoms_b->GetPosition()));
          if( sqr_distance < 0.1 || sqr_distance > dist_cutoff_sqr)
          {
            // overlapping atom
            continue;
          }
          force_sum += transition( sqr_distance) * *itr_charge_b * *itr_charge_a / sqr_distance;
        }
      }
      STORAGE( 0) = force_sum;
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CoulombicForce::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes the intramolecular coulombic force (without the dielectric constant)" +
        std::string(
          m_ExcludeNeighborsOfNeighbors
          ? " due to the rotamer chosen, thus ignoring covalent and 1,3 bonding interactions"
          : ""
        )
      );

      parameters.AddInitializer
      (
        "",
        "charge property",
        io::Serialization::GetAgent( &m_Charge)
      );
      parameters.AddInitializer
      (
        "distance cutoff",
        "maximum distance to consider. forces will be scaled (using a cosine transition) from 1A prior to this value",
        io::Serialization::GetAgent( &m_DistanceCutoff),
        "nan"
      );
      if( !m_ExcludeNeighborsOfNeighbors)
      {
        parameters.AddInitializer
        (
          "non-covalent",
          "true to ignore covalently bonded neighbors of an atom for the force calculation",
          io::Serialization::GetAgent( &m_ExcludeNeighbors),
          "False"
        );
      }
      return parameters;
    }
  } // namespace descriptor
} // namespace bcl
