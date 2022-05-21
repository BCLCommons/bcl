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
#include "descriptor/bcl_descriptor_atom_aromaticity_axes.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_conformational.h"
#include "chemistry/bcl_chemistry_constitutional_bond_type_data.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    AtomAromaticityAxes *AtomAromaticityAxes::Clone() const
    {
      return new AtomAromaticityAxes( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomAromaticityAxes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomAromaticityAxes::GetAlias() const
    {
      static const std::string s_Name( "Atom_AromaticityAxes");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomAromaticityAxes::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      size_t n_aro_neigh( ELEMENT->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsAromatic, 1));
      if( !n_aro_neigh)
      {
        return;
      }
      storage::Vector< linal::Vector3D> aro_points( n_aro_neigh + 1);
      aro_points( 0) = ELEMENT->GetPosition();
      size_t pos( 1);
      for( auto itr( ELEMENT->GetBonds().Begin()), itr_end( ELEMENT->GetBonds().End()); itr != itr_end; ++itr)
      {
        if( itr->GetBondType()->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Aromatic)
        {
          aro_points( pos++) = itr->GetTargetAtom().GetPosition();
        }
      }
      linal::Vector3D ab( aro_points( 1) - aro_points( 0)), bc( aro_points( 2) - aro_points( 0));
      linal::Vector3D cross_product( linal::CrossProduct( ab, bc));
      cross_product.Normalize();
      double cross_sum
      (
        ( cross_product.X() < 0.0 ? -1 : 1) * math::Sqr( cross_product.X())
        + ( cross_product.Y() < 0.0 ? -1 : 1) * math::Sqr( cross_product.Y())
        + ( cross_product.Z() < 0.0 ? -1 : 1) * math::Sqr( cross_product.Z())
      );
      if( cross_sum < 0.0)
      {
        cross_product = -cross_product;
      }
      std::copy( cross_product.Begin(), cross_product.End(), STORAGE.Begin());
    } // Recalculate

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomAromaticityAxes::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "For aromatic atoms, the axes of the aromatic field");
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
