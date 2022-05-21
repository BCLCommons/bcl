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
#include "descriptor/bcl_descriptor_atom_polarizability.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    AtomPolarizability *AtomPolarizability::Clone() const
    {
      return new AtomPolarizability( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomPolarizability::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomPolarizability::GetAlias() const
    {
      static const std::string s_Name( "Atom_Polarizability");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the polarizability for a particular atom
    //! @param ATOM the atom of interest
    //! @return the polarizability of the given atom
    float AtomPolarizability::GetPolarizability( const chemistry::AtomConformationalInterface &ATOM)
    {
      // get the basic additive polarizability for the atom type
      float polarizability
      (
        ATOM.GetAtomType()->GetAtomTypeProperty( chemistry::AtomTypeData::e_AdditiveAtomicPolarizability)
      );

      // for C_TrTrTrPi, change the value to 1.896 if not connected to any hydrogens
      //                 otherwise, set it to 1.352
      // see J.Am.Chem.Soc. Vol 112, No. 23, 1990, 8534
      if( ATOM.GetAtomType() == chemistry::GetAtomTypes().C_TrTrTrPi)
      {
        // initialize with whether there are implicit hydrogens
        bool is_bonded_to_h
        (
          ATOM.GetNumberCovalentlyBoundHydrogens() ||
          ATOM.CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_BondOrder, 1) != size_t( 2)
        );
        polarizability = is_bonded_to_h ? 1.352 : 1.896;
      }
      else if( !util::IsDefined( polarizability))
      {
        // undefined polarizability; do what adriana does and just use the value for TeTeTeTe
        polarizability = chemistry::GetAtomTypes().C_TeTeTeTe->GetAtomTypeProperty( chemistry::AtomTypeData::e_AdditiveAtomicPolarizability);
      }
      return polarizability;
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomPolarizability::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      STORAGE( 0) = GetPolarizability( *ELEMENT);
    } // Recalculate

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomPolarizability::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes the polarizability of each atom using the method from see J.Am.Chem.Soc. Vol 112, No. 23, 1990, 8534"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
