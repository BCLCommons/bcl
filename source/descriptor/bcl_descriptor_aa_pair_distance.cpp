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
#include "descriptor/bcl_descriptor_aa_pair_distance.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAPairDistance::s_Instance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AAPairDistance()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AAPairDistance *AAPairDistance::Clone() const
    {
      return new AAPairDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAPairDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAPairDistance::GetAlias() const
    {
      static const std::string s_name( "AAPairDistance");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAPairDistance::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT_A, ELEMENT_B: the element pair of interest
    //! @param STORAGE storage for the descriptor
    void AAPairDistance::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT_A,
      const iterate::Generic< const biol::AABase> &ELEMENT_B,
      linal::VectorReference< float> &STORAGE
    )
    {
      STORAGE( 0) = biol::FirstSidechainAtomDistance( *ELEMENT_A, *ELEMENT_B);
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPairDistance::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Returns the distance between the amino acid pair in Angstrom");

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
