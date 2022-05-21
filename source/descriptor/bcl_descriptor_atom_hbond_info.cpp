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
#include "descriptor/bcl_descriptor_atom_hbond_info.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomHBondInfo::AtomHBondInfo( const Method &METHOD) :
        m_Type( METHOD)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new AtomHBondInfo
    AtomHBondInfo *AtomHBondInfo::Clone() const
    {
      return new AtomHBondInfo( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomHBondInfo::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomHBondInfo::GetAlias() const
    {
      static std::string s_hbond_acc( "Atom_HbondAcceptors");
      static std::string s_hbond_don( "Atom_HbondDonors");
      return m_Type == e_Acceptor ? s_hbond_acc : s_hbond_don;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomHBondInfo::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // do we have an oxygen or a nitrogen?
      STORAGE( 0) = 0.0;
      if
      (
        ELEMENT->GetElementType() == chemistry::GetElementTypes().e_Oxygen
        || ELEMENT->GetElementType() == chemistry::GetElementTypes().e_Nitrogen
      )
      {
        if( m_Type == e_Acceptor)
        {
          // yep, so the atom is a hydrogen bond acceptor
          STORAGE( 0) = 1.0;
        }
        // check if any of the atoms are hydrogen
        else if( ELEMENT->GetNumberCovalentlyBoundHydrogens() > size_t( 0))
        {
          // yep, so the atom is a hydrogen bond donor
          STORAGE( 0) = 1.0;
        }
        // if there are any single bond valences, then count this atom as well as an hbond donor
        else if( ELEMENT->GetNumberElectronsInValenceBonds() < 2 * ELEMENT->GetNumberValenceBonds())
        {
          STORAGE( 0) = 1.0;
        }
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomHBondInfo::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_Type == e_Acceptor ? "1 for hydrogen bond acceptors (N and O), 0 for other elements" :
            "1 for hydrogen bond donors (NH and OH), 0 for others"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
