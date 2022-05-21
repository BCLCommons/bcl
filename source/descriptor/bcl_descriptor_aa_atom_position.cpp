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
#include "descriptor/bcl_descriptor_aa_atom_position.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AAAtomPosition::s_AATMPositionInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AAAtomPosition()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAAtomPosition::AAAtomPosition() :
      m_Type( biol::GetAtomTypes().CA)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AAAtomPosition *AAAtomPosition::Clone() const
    {
      return new AAAtomPosition( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAAtomPosition::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAAtomPosition::GetAlias() const
    {
      static const std::string s_name( "AA_AtomPosition");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAAtomPosition::GetNormalSizeOfFeatures() const
    {
      return 3;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT_A the element of interest
    //! @param STORAGE storage for the descriptor
    void AAAtomPosition::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      const linal::Vector3D &position( ELEMENT->GetAtom( m_Type).GetCoordinates());
      std::copy( position.Begin(), position.End(), STORAGE.Begin());
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAAtomPosition::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Retrieves the position of the given atom in the aa");

      parameters.AddInitializer
      (
        "",
        "atom type for which the position should be retrieved",
        io::Serialization::GetAgent( &m_Type),
        "CA"
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
