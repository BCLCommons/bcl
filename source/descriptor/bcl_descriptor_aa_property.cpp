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
#include "descriptor/bcl_descriptor_aa_property.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for( size_t aa_property( 0); aa_property < biol::AATypeData::s_NumberPropertyTypes; ++aa_property)
        {
          last_instance =
            util::Enumerated< Base< biol::AABase, float> >::AddInstance
            (
              new AAProperty( static_cast< biol::AATypeData::PropertyType>( aa_property))
            );
        }
        return last_instance;
      }
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAProperty::s_Instance( AddInstances());

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a list of properties
    //! @param PROPERTY property to be used
    AAProperty::AAProperty
    (
      const biol::AATypeData::PropertyType &PROPERTY
    ) :
      m_Property( PROPERTY)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AABlastProfile
    AAProperty *AAProperty::Clone() const
    {
      return new AAProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAProperty::GetAlias() const
    {
      return m_Property.GetString();
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAProperty::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AAProperty::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // if we have values for the aa properties, use them, otherwise, use the values from the arbitrary amino acid
      if
      (
        ELEMENT->GetType() == biol::GetAATypes().R1A
        || ELEMENT->GetType() == biol::GetAATypes().UNK
        || ELEMENT->GetType() == biol::GetAATypes().GAP
        || !ELEMENT->GetType().IsDefined()
      )
      {
        STORAGE( 0) = biol::GetAATypes().XXX->GetAAProperty( m_Property);
      }
      else
      {
        STORAGE( 0) = ELEMENT->GetType()->GetAAProperty( m_Property);

        // TPSA varies substantially at the terminii, unlike most other properties, so modify its return value
        // if at either terminus
        if( m_Property == biol::AATypeData::e_SideChainTopologicalPolarSurfaceArea)
        {
          if( ELEMENT.GetPosition() == 0)
          {
            // N-terminus, add 27.64, except for proline (which has a different N atom type, hence a lower polarizability)
            if( ELEMENT->GetType() == biol::GetAATypes().PRO)
            {
              STORAGE( 0) += 16.61;
            }
            else
            {
              STORAGE( 0) += 27.64;
            }
          }
          if( ELEMENT.GetReversePosition() == 1)
          {
            // C-terminus, add 23 for terminal carboxylate
            STORAGE( 0) += 23.0;
          }
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAProperty::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( GetAlias() + " for each amino acid");
      return serializer;
    }

  } // namespace descriptor
} // namespace bcl
