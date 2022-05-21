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
#include "descriptor/bcl_descriptor_molecule_girth.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeGirth
    MoleculeGirth *MoleculeGirth::Clone() const
    {
      return new MoleculeGirth( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeGirth::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeGirth::GetAlias() const
    {
      static const std::string s_name( "Girth");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeGirth::Calculate( linal::VectorReference< float> &STORAGE)
    {
      float max_distance( 0);
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_a( this->GetCurrentObject()->GetIterator());
        itr_atoms_a.NotAtEnd();
        ++itr_atoms_a
      )
      {
        for
        (
          iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_b( itr_atoms_a);
          itr_atoms_b.NotAtEnd();
          ++itr_atoms_b
        )
        {
          // store distance between both atoms
          const float distance( linal::Distance( itr_atoms_a->GetPosition(), itr_atoms_b->GetPosition()));

          // update max distance
          if( distance > max_distance)
          {
            max_distance = distance;
          }
        }
      }
      STORAGE = max_distance;
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeGirth::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "calculates the girth of a molecule");
      return parameters;
    }
  } // namespace descriptor
} // namespace bcl
