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
#include "descriptor/bcl_descriptor_aa_identifier.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "pdb/bcl_pdb_entry_types.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAIdentifier::s_SeqIDInstance
    (
      util::Enumerated< Base< biol::AABase, char> >::AddInstance
      (
        new AAIdentifier( true)
      )
    );
    const util::SiPtr< const util::ObjectInterface> AAIdentifier::s_PdbIDInstance
    (
      util::Enumerated< Base< biol::AABase, char> >::AddInstance
      (
        new AAIdentifier( false)
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor with ID type
    AAIdentifier::AAIdentifier( const bool &SEQ_ID) :
      m_SeqId( SEQ_ID)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AAIdentifier *AAIdentifier::Clone() const
    {
      return new AAIdentifier( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAIdentifier::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAIdentifier::GetAlias() const
    {
      static const std::string s_seq_name( "AASeqID"), s_pdb_name( "AAPdbID");
      return m_SeqId ? s_seq_name : s_pdb_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAIdentifier::GetNormalSizeOfFeatures() const
    {
      return pdb::GetEntryTypes().ATOMSegmentID->GetLength() + 1;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AAIdentifier::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< char> &STORAGE
    )
    {
      const std::string id
      (
        util::Format().W( GetNormalSizeOfFeatures() - 1)( m_SeqId ? ELEMENT->GetSeqID() : ELEMENT->GetPdbID())
      );
      // always prepend with a space for readability
      STORAGE( 0) = ' ';
      std::copy( id.begin(), id.end(), STORAGE.Begin() + 1);
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAIdentifier::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "the amino acid's " + std::string( m_SeqId ? "sequence" : "PDB") + " ID");

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
