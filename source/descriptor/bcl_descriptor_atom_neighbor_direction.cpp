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
#include "descriptor/bcl_descriptor_atom_neighbor_direction.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_conformational.h"
#include "chemistry/bcl_chemistry_constitutional_bond_type_data.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_running_average.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    AtomNeighborDirection *AtomNeighborDirection::Clone() const
    {
      return new AtomNeighborDirection( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomNeighborDirection::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomNeighborDirection::GetAlias() const
    {
      static const std::string s_Name( "Atom_NeighborDirection");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomNeighborDirection::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      //Calculate average position unit vector in plane of neighboring atoms
      math::RunningAverage< linal::Vector3D> neighbor_avg;
      for( auto itr( ELEMENT->GetBonds().Begin()), itr_end( ELEMENT->GetBonds().End()); itr < itr_end; ++itr)
      {
        neighbor_avg += linal::UnitVector( ELEMENT->GetPosition(), itr->GetTargetAtom().GetPosition());
      }
      std::copy( neighbor_avg.GetAverage().Begin(), neighbor_avg.GetAverage().End(), STORAGE.Begin());
      float norm( STORAGE.Norm());
      if( norm > 1e-3)
      {
        STORAGE /= norm;
      }
    }
    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomNeighborDirection::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "The unit vector formed by the relative positions of neighboring atoms");
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
