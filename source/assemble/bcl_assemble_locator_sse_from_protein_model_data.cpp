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
#include "assemble/bcl_assemble_locator_sse_from_protein_model_data.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorSSEFromProteinModelData::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorSSEFromProteinModelData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorSSEFromProteinModelData::LocatorSSEFromProteinModelData() :
      m_Key( ProteinModelData::e_Undefined)
    {
    }

    //! @brief constructor from a key
    //! @param KEY key to retrieve the locators from the protein model data
    LocatorSSEFromProteinModelData::LocatorSSEFromProteinModelData( const ProteinModelData::Type KEY) :
      m_Key( KEY)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LocatorSSEFromProteinModelData
    LocatorSSEFromProteinModelData *LocatorSSEFromProteinModelData::Clone() const
    {
      return new LocatorSSEFromProteinModelData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorSSEFromProteinModelData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locate and return an SSE using the locators retrieved from ProteinModelData
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return returns SiPtr to selected SSE
    util::SiPtr< const SSE> LocatorSSEFromProteinModelData::Locate( const DomainInterface &PROTEIN_MODEL) const
    {
      // try to cast the DomainInterface to ProteinModel
      util::SiPtr< const ProteinModel> sp_model( &PROTEIN_MODEL);
      if( !sp_model.IsDefined())
      {
        BCL_MessageCrt( "The cast from DomainInterface to ProteinModel failed!");
        return util::SiPtr< const SSE>();
      }

      // get pointer to data from protein model
      const util::ShPtr< util::ShPtrList< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > > sp_data
      (
        sp_model->GetProteinModelData()->GetData( m_Key)
      );

      // if not defined then warn user and return empty SiPtr
      if( !sp_data.IsDefined())
      {
        BCL_MessageCrt
        (
          "Could not locate protein model data with key " + ProteinModelData::GetTypeName( m_Key)
        );
        return util::SiPtr< const SSE>();
      }

      // if the size is equal to 0 meaning there are no such locators
      if( sp_data->IsEmpty())
      {
        return util::SiPtr< const SSE>();
      }

      // get a random iterator
      util::ShPtrList< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> >::const_iterator itr
      (
        random::GetGlobalRandom().Iterator( sp_data->Begin(), sp_data->End(), sp_data->GetSize())
      );

      // use the randomly located locator to locate an SSE and return it
      return ( *itr)->Locate( PROTEIN_MODEL);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorSSEFromProteinModelData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Key, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorSSEFromProteinModelData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Key, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
