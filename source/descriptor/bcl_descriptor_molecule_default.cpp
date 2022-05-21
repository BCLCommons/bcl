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
#include "descriptor/bcl_descriptor_molecule_default.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeDefault::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeDefault( "Default")
      )
    );

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeDefault::s_DefaultInstance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeDefault( "")
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeDefault::MoleculeDefault() :
      m_Alias( ""),
      m_MiscPropertyString(),
      m_Constant( util::GetUndefined< float>())
    {
    }

    //! @brief constructor from a alias
    //! @param ALIAS the alias to use
    MoleculeDefault::MoleculeDefault( const std::string &ALIAS) :
      m_Alias( ALIAS),
      m_MiscPropertyString(),
      m_Constant()
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeDefault
    MoleculeDefault *MoleculeDefault::Clone() const
    {
      return new MoleculeDefault( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeDefault::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeDefault::GetAlias() const
    {
      return m_Alias;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeDefault::Calculate( linal::VectorReference< float> &STORAGE)
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);

      // check for whether a numerical value was given
      if( util::IsDefined( m_Constant))
      {
        STORAGE = m_Constant;
        return;
      }

      // check for the property in the cache
      if( molecule.IsCached( m_MiscPropertyString))
      {
        STORAGE.CopyValues( molecule.GetFromCache( m_MiscPropertyString));
      }
      // check for the property in the storage
      else if( molecule.IsPropertyStored( m_MiscPropertyString))
      {
        STORAGE.CopyValues( molecule.GetMDLPropertyAsVector( m_MiscPropertyString));
      }
      else
      {
        STORAGE = util::GetUndefined< float>();
        BCL_MessageStd( " no values found from default " + util::Format()( m_MiscPropertyString));
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeDefault::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( m_MiscPropertyString.empty())
      {
        if( LABEL.IsScalar() && !LABEL.GetValue().empty())
        {
          m_MiscPropertyString = LABEL.GetValue();
        }
        else
        {
          ERR_STREAM << "No property or number was given!";
          return false;
        }
      }

      if( util::IsNumerical( m_MiscPropertyString))
      {
        m_Constant = util::ConvertStringToNumericalValue< float>( m_MiscPropertyString);
      }
      else
      {
        m_Constant = util::GetUndefined< float>();
      }
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeDefault::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "given a number, returns that number, otherwise retrieves the value of the misc property by that name"
      );

      parameters.AddInitializer
      (
        "",
        "number or misc property name",
        io::Serialization::GetAgent( &m_MiscPropertyString),
        ""
      );

      return parameters;
    }
  } // namespace descriptor
} // namespace bcl
