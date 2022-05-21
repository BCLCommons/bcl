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
#include "descriptor/bcl_descriptor_atom_volume.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_conformational.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomVolume::AtomVolume()
    {
    }

    //! @brief default constructor
    //! @param COVALENT_RADIUS whether to use covalent radius
    //! @param USE_CSD_VDW only used if
    AtomVolume::AtomVolume
    (
      const util::Implementation< Base< chemistry::AtomConformationalInterface, float> > &MAX_RADIUS,
      const util::Implementation< Base< chemistry::AtomConformationalInterface, float> > &MIN_RADIUS
    ) :
      m_MaxRadius( MAX_RADIUS),
      m_MinRadius( MIN_RADIUS)
    {
    }

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AtomVolume::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new AtomVolume()
      )
    );

    //! @brief virtual copy constructor
    //! @return pointer to new AtomVolume
    AtomVolume *AtomVolume::Clone() const
    {
      return new AtomVolume( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomVolume::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomVolume::GetAlias() const
    {
      static std::string s_name( "Atom_Volume");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > AtomVolume::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_MaxRadius,
          m_MinRadius.IsDefined() ? &m_MinRadius + 1 : &m_MaxRadius + 1
        );
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomVolume::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      Iterator< chemistry::AtomConformationalInterface> itr( ELEMENT);
      float max_radii( m_MaxRadius->operator ()( itr)( 0));
      float min_radii( m_MinRadius.IsDefined() ? m_MinRadius->operator ()( itr)( 0) : 0.0);

      // start off with the volume == to the max volume
      float volume( 4.0 / 3.0 * math::g_Pi * math::Pow( max_radii, float( 3.0)));

      // each atom has at least the min volume
      const float min_volume( 4.0 / 3.0 * math::g_Pi * math::Pow( min_radii, float( 3.0)));

      for
      (
        storage::Vector< chemistry::BondConformational>::const_iterator
          itr_conn( ELEMENT->GetBonds().Begin()), itr_conn_end( ELEMENT->GetBonds().End());
        itr_conn != itr_conn_end && volume > min_volume;
        ++itr_conn
      )
      {
        // calculate the distance between the current atom and current bonded atom.
        // to calculate distance based on 3d coordinates of atoms.
        const float distance
        (
          linal::Distance( ELEMENT->GetPosition(), itr_conn->GetTargetAtom().GetPosition())
        );

        Iterator< chemistry::AtomConformationalInterface> itr_bonded
        (
          iterate::Generic< const chemistry::AtomConformationalInterface>
          (
            &itr_conn->GetTargetAtom(),
            &itr_conn->GetTargetAtom() + 1
          )
        );
        float max_radius_bonded_atom( m_MaxRadius->operator ()( itr_bonded)( 0));

        // detect distances which indicate that the atoms overlap completely or not at all
        if( distance == 0.0 || distance >= max_radii + max_radius_bonded_atom)
        {
          continue;
        }

        // Remove from the current atom's surface area the calculated overlap between
        // the current atom and the current bonded atom.
        // This uses a simple spherical overlap equation:
        // radius1 * pi * [(radius2^2 - (radius1 - distance)^2) / distance]
        // Formula can be found in J. Mol Graphics and Modeling 18. Author Paul Labute.
        // this calculation assumes that no more than two atoms
        const float volume_cap_lost
        (
          math::g_Pi * math::Sqr( max_radii + max_radius_bonded_atom - distance)
          *
          (
            math::Sqr( distance) + 2 * distance * ( max_radius_bonded_atom + max_radii)
            - 3 * ( math::Sqr( max_radius_bonded_atom) + math::Sqr( max_radii)) + 6 * max_radius_bonded_atom * max_radii
          ) / ( 12 * distance)
        );
        volume -= volume_cap_lost;
      }
      // Once all of the bonded atoms have been iterated, add the remaining surface area to the vector
      // of all surface areas in the molecule.
      STORAGE( 0) = std::max( volume, min_volume);

    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomVolume::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "approximates the volume of the atom, considering neighbor overlap");
      parameters.AddInitializer
      (
        "radius",
        "Descriptor that defines the maximum atomic radius, assuming no overlap from neighboring atoms",
        io::Serialization::GetAgent( &m_MaxRadius)
      );
      parameters.AddInitializer
      (
        "min radius",
        "Descriptor that defines the minimum atomic radius, after accounting for overlap from neighboring atoms",
        io::Serialization::GetAgent( &m_MinRadius),
        "0"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
