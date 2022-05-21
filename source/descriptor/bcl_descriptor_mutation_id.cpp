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
#include "descriptor/bcl_descriptor_mutation_id.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
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
    const util::SiPtr< const util::ObjectInterface> MutationId::s_Instance
    (
      util::Enumerated< Base< biol::Mutation, char> >::AddInstance( new MutationId())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new BaseElement
    MutationId *MutationId::Clone() const
    {
      return new MutationId( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutationId::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &MutationId::GetAlias() const
    {
      static const std::string s_name( "MutationId");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MutationId::GetNormalSizeOfFeatures() const
    {
      return 6;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MutationId::Calculate
    (
      const iterate::Generic< const biol::Mutation> &ELEMENT,
      linal::VectorReference< char> &STORAGE
    )
    {
      // copy the cached protein id
      std::string id( ELEMENT->ToString());
      while( id.size() < GetNormalSizeOfFeatures())
      {
        id += ' ';
      }
      STORAGE.CopyValues
      (
        linal::VectorConstReference< char>( GetNormalSizeOfFeatures(), &id[ 0])
      );
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutationId::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "the mutation id, written as wild-type amino acid (1-letter code), SeqID, mutant aa type (1-letter code)");

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
