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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "descriptor/bcl_descriptor_cache_map.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_cache_map.cpp
  //!
  //! @author mendenjl
  //! @date Nov 06, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorCacheMap :
    public ExampleInterface
  {
  public:

    ExampleDescriptorCacheMap *Clone() const
    {
      return new ExampleDescriptorCacheMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      descriptor::CacheMap props;

      // test that does mdl property exists does not return false positives
      BCL_ExampleCheck
      (
        descriptor::CacheMap().Has( util::ObjectDataLabel( "property")),
        false
      );

      // test that does mdl property exists does not return false negatives
      props.Insert( util::ObjectDataLabel( "property"), linal::Vector< float>( 1, float( 1.0)));
      BCL_ExampleIndirectAssert
      (
        props.Has( util::ObjectDataLabel( "property")),
        true,
        "StoreProperty"
      );

      // test that we get the property string back properly
      BCL_ExampleCheck
      (
        props.Get( util::ObjectDataLabel( "property")),
        linal::Vector< float>( 1, float( 1.0))
      );

      props.Insert( util::ObjectDataLabel( "property"), linal::Vector< float>( 1, float( 2.0)));

      // check that there is really just one mdl property in the map
      BCL_ExampleIndirectCheck( props.GetSize(), 1, "Entries in CacheMap should be overwriteable");

      // check that we can add a property
      props.Insert( util::ObjectDataLabel( "another property"), storage::Vector< float>( 1, float( 3.5)));
      BCL_ExampleIndirectCheck
      (
        props.GetSize(),
        2,
        "CacheMap with second unique property name"
      );

      // check that the new property and the original property have the correct values
      BCL_ExampleIndirectCheck
      (
        props.Get( util::ObjectDataLabel( "property")),
        linal::Vector< float>( 1, float( 2.0)),
        "CacheMap(a) should not affect property b"
      );
      BCL_ExampleIndirectCheck
      (
        props.Get( util::ObjectDataLabel( "another property")),
        linal::Vector< float>( 1, float( 3.5)),
        "CacheMap"
      );

      props.Reset();
      BCL_ExampleIndirectCheck( props.GetSize(), 0, "Reset");

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  };

  const ExampleClass::EnumType ExampleDescriptorCacheMap::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorCacheMap())
  );

} // namespace bcl
