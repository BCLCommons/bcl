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
#include "descriptor/bcl_descriptor_protein_id.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ProteinId::s_Instance
    (
      util::Enumerated< Base< biol::AABase, char> >::AddInstance( new ProteinId())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new BaseElement
    ProteinId *ProteinId::Clone() const
    {
      return new ProteinId( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ProteinId::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &ProteinId::GetAlias() const
    {
      static const std::string s_name( "ProteinId");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t ProteinId::GetNormalSizeOfFeatures() const
    {
      return m_NumChar;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void ProteinId::Calculate( linal::VectorReference< char> &STORAGE)
    {
      // copy the cached protein id
      STORAGE.CopyValues( linal::VectorConstReference< char>( GetNormalSizeOfFeatures(), &m_Id[ 0]));
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinId::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "the PDB ID of the protein");
      parameters.AddOptionalInitializer
      (
        "",
        "Number of characters to include for each ID. For PDB-id'ed files, the default (5) is appropriate for including "
        "chain ids",
        io::Serialization::GetAgent( &m_NumChar)
      );

      return parameters;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void ProteinId::SetObjectHook()
    {
      util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model( this->GetCurrentObject());

      // Get the filename
      util::ShPtr< util::Wrapper< std::string> > sp_id_wrapper
      (
        sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_Identification)
      );

      // get a reference to the string
      const std::string &id_str( sp_id_wrapper->GetData());

      // handle excessively long ids
      if( id_str.size() > GetNormalSizeOfFeatures())
      {
        m_Id = id_str.substr( 0, GetNormalSizeOfFeatures());
        BCL_MessageVrb( "Excessively long id: " + id_str + " pruned to " + m_Id);
      }
      else
      {
        m_Id = id_str;
        m_Id.resize( GetNormalSizeOfFeatures(), ' ');
      }
    }

  } // namespace descriptor
} // namespace bcl
