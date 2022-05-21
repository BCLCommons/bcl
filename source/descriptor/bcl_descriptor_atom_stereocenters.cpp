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
#include "descriptor/bcl_descriptor_atom_stereocenters.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_stereocenters_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new AtomStereocenters
    AtomStereocenters *AtomStereocenters::Clone() const
    {
      return new AtomStereocenters( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomStereocenters::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomStereocenters::GetAlias() const
    {
      static const std::string s_name( "Atom_Stereocenters");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomStereocenters::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      if( m_Stereocenters.IsEmpty())
      {
        RecalculateStereocenters();
      }
      STORAGE( 0) = m_Stereocenters( ELEMENT.GetPosition());
    } // Recalculate

    //! @brief Recalculates the charges for the current molecule
    void AtomStereocenters::RecalculateStereocenters()
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);

      m_Stereocenters = chemistry::StereocentersHandler::CalculateFromConformation( molecule.GetAtomsIterator());
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AtomStereocenters::SetObjectHook()
    {
      m_Stereocenters = linal::Vector< float>();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomStereocenters::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "1 for R, -1 for S, 0 for achiral atoms, 2 for undefined chirality");
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
