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
#include "descriptor/bcl_descriptor_atom_estimated_surface_area.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param COVALENT_RADIUS whether to use covalent radius
    //! @param USE_CSD_VDW only used if
    AtomEstimatedSurfaceArea::AtomEstimatedSurfaceArea
    (
      const bool &COVALENT_RADIUS,
      const bool &USE_CSD_VDW
    ) :
      m_UseCovalentRadius( COVALENT_RADIUS),
      m_UseCSDVdwRadius( USE_CSD_VDW)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new AtomEstimatedSurfaceArea
    AtomEstimatedSurfaceArea *AtomEstimatedSurfaceArea::Clone() const
    {
      return new AtomEstimatedSurfaceArea( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomEstimatedSurfaceArea::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomEstimatedSurfaceArea::GetAlias() const
    {
      static std::string s_covalent_name( "Atom_EstCovalentSurfaceArea"), s_vdw_name( "Atom_EstVdwSurfaceArea");
      static std::string  s_csd_vdw_name( "Atom_EstVdwSurfaceAreaCSD");
      return m_UseCSDVdwRadius && !m_UseCovalentRadius ? s_csd_vdw_name : ( m_UseCovalentRadius ? s_covalent_name : s_vdw_name);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomEstimatedSurfaceArea::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {

      // calculate the covalent surface area
      const float ave_covalent_radius( chemistry::BondLengths::GetAverageCovalentRadius( *ELEMENT));

      // Calculate the surface area of the current atom. 4*pi* radius^2
      const float covalent_surface_area( 4 * math::g_Pi * math::Sqr( ave_covalent_radius));

      // if just calculating covalent radii, this is all that we have to do
      if( m_UseCovalentRadius)
      {
        STORAGE( 0) = covalent_surface_area;
        return;
      }

      // van der waals radius; use sphere overlap equation

      // calculate the van-der-waals surface area of the atom
      const float vdw_radius
      (
        m_UseCSDVdwRadius
        ? ELEMENT->GetAtomType()->GetAtomTypeProperty( chemistry::AtomTypeData::e_VdWaalsRadiusCSD)
        : ELEMENT->GetElementType()->GetProperty( chemistry::ElementTypeData::e_VDWaalsRadius)
      );
      float vdw_surface_area( 4 * math::g_Pi * math::Sqr( vdw_radius));
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
        const float covalent_distance
        (
          chemistry::BondLengths::GetBondLength
          (
            ELEMENT->GetAtomType(),
            itr_conn->GetBondType()->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromatic),
            itr_conn->GetTargetAtom().GetAtomType()
          )
        );
        const float vdw_radius_connection
        (
          m_UseCSDVdwRadius
          ? itr_conn->GetTargetAtom().GetAtomType()->GetAtomTypeProperty( chemistry::AtomTypeData::e_VdWaalsRadiusCSD)
          : itr_conn->GetTargetAtom().GetElementType()->GetProperty( chemistry::ElementTypeData::e_VDWaalsRadius)
        );

        // detect distances which indicate that the atoms overlap completely or not at all
        if
        (
          covalent_distance <= math::Absolute( vdw_radius - vdw_radius_connection)
          || covalent_distance >= vdw_radius + vdw_radius_connection
        )
        {
          continue;
        }

        // Remove from the current atom's surface area the calculated overlap between
        // the current atom and the current bonded atom.
        // This uses a simple spherical overlap equation:
        // radius1 * pi * [(radius2^2 - (radius1 - distance)^2) / distance]
        // Formula can be found in J. Mol Graphics and Modeling 18. Author Paul Labute.
        // this calculation assumes that no more than two atoms
        vdw_surface_area -= math::g_Pi * vdw_radius * ( math::Sqr( vdw_radius_connection) - math::Sqr( vdw_radius - covalent_distance)) / covalent_distance;
      }
      // Once all of the bonded atoms have been iterated, add the remaining surface area to the vector
      // of all surface areas in the molecule.
      STORAGE( 0) = std::max( covalent_surface_area, vdw_surface_area);

    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomEstimatedSurfaceArea::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "The surface area, estimated from atom type "
        + std::string
          (
            m_UseCovalentRadius
            ? "covalent"
            : ( m_UseCSDVdwRadius ? "CSD-derived van der waals" : "elemental van der waals")
          )
        + " radii"
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
