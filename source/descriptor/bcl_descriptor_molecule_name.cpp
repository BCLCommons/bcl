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
#include "descriptor/bcl_descriptor_molecule_name.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // single instance of that class
    const util::SiPtr< const util::ObjectInterface> MoleculeName::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, char> >::AddInstance
      (
        new MoleculeName()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeName::MoleculeName() :
      m_NumberCharacters( 240)
    {
    }

    //! virtual copy constructor
    MoleculeName *MoleculeName::Clone() const
    {
      return new MoleculeName( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoleculeName::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoleculeName::GetAlias() const
    {
      static const std::string s_name( "Name");
      return s_name;
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    void MoleculeName::RecalculateImpl
    (
      const Iterator< chemistry::AtomConformationalInterface> &ITR,
      linal::VectorReference< char> &STORAGE
    )
    {
      util::SiPtr< const chemistry::ConformationInterface> molecule_pointer( GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *molecule_pointer);
      std::string original_name( util::TrimString( molecule.GetName()));
      std::replace( original_name.begin(), original_name.end(), '\n', '\t');
      // return the name, with newlines replaced by tabs

      if( original_name.size() <= m_NumberCharacters)
      {
        std::copy( original_name.begin(), original_name.end(), STORAGE.Begin());
      }
      else
      {
        std::copy( original_name.begin(), original_name.begin() + m_NumberCharacters, STORAGE.Begin());
      }

    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeName::GetSerializer() const
    {
      io::Serializer init;
      init.SetClassDescription( "Retrieves the name of the molecule");
      init.AddInitializer
      (
        "size",
        "number of characters expected.  If the string is longer than this, it will be truncated to this length."
        "If it is shorter, it will be padded with spaces to this number of characters, provided that it exists on the "
        "molecule (otherwise it will be blank)",
        io::Serialization::GetAgent( &m_NumberCharacters),
        "240"
      );

      return init;
    }
  } // namespace descriptor
} // namespace bcl
