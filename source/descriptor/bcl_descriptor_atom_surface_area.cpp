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
#include "descriptor/bcl_descriptor_atom_surface_area.h"

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
    AtomSurfaceArea::AtomSurfaceArea()
    {
    }

    //! @brief default constructor
    //! @param COVALENT_RADIUS whether to use covalent radius
    //! @param USE_CSD_VDW only used if
    AtomSurfaceArea::AtomSurfaceArea
    (
      const util::Implementation< Base< chemistry::AtomConformationalInterface, float> > &MAX_RADIUS,
      const util::Implementation< Base< chemistry::AtomConformationalInterface, float> > &MIN_RADIUS
    ) :
      m_MaxRadius( MAX_RADIUS),
      m_MinRadius( MIN_RADIUS)
    {
    }

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AtomSurfaceArea::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new AtomSurfaceArea()
      )
    );

    //! @brief virtual copy constructor
    //! @return pointer to new AtomSurfaceArea
    AtomSurfaceArea *AtomSurfaceArea::Clone() const
    {
      return new AtomSurfaceArea( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomSurfaceArea::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomSurfaceArea::GetAlias() const
    {
      static std::string s_name( "Atom_SurfaceArea");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > AtomSurfaceArea::GetInternalDescriptors()
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
    void AtomSurfaceArea::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      Iterator< chemistry::AtomConformationalInterface> itr( ELEMENT);
      float max_radii( m_MaxRadius->operator ()( itr)( 0));
      float min_radii( m_MinRadius.IsDefined() ? m_MinRadius->operator ()( itr)( 0) : 0.0);

      // Calculate the surface area of the current atom. 4*pi* radius^2
      float surface_area( 4 * math::g_Pi * math::Sqr( max_radii));

      // Calculate the minimum surface area for the current atom
      const float minimum_surface_area( 4 * math::g_Pi * math::Sqr( min_radii));

      for
      (
        storage::Vector< chemistry::BondConformational>::const_iterator
          itr_conn( ELEMENT->GetBonds().Begin()), itr_conn_end( ELEMENT->GetBonds().End());
        itr_conn != itr_conn_end;
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
        if( distance <= 0.0 || distance >= max_radii + max_radius_bonded_atom)
        {
          continue;
        }

        // Remove from the current atom's surface area the calculated overlap between
        // the current atom and the current bonded atom.
        // This uses a simple spherical overlap equation:
        // radius1 * pi * [(radius2^2 - (radius1 - distance)^2) / distance]
        // Formula can be found in J. Mol Graphics and Modeling 18. Author Paul Labute.
        // this calculation assumes that no more than two atoms
        surface_area -= math::g_Pi * max_radii * ( math::Sqr( max_radius_bonded_atom) - math::Sqr( max_radii - distance)) / distance;
      }
      // Once all of the bonded atoms have been iterated, add the remaining surface area to the vector
      // of all surface areas in the molecule.
      STORAGE( 0) = std::max( surface_area, minimum_surface_area);

    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomSurfaceArea::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "approximates the surface area of the atom, considering neighbor overlap");

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
