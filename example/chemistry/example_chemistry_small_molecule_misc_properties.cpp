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
#include "chemistry/bcl_chemistry_small_molecule_misc_properties.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_small_molecule_misc_properties.cpp
  //!
  //! @author mendenjl
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistrySmallMoleculeMiscProperties :
    public ExampleInterface
  {
  public:

    ExampleChemistrySmallMoleculeMiscProperties *Clone() const
    {
      return new ExampleChemistrySmallMoleculeMiscProperties( *this);
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

      chemistry::SmallMoleculeMiscProperties props;

      // test that does mdl property exists does not return false positives
      BCL_ExampleCheck
      (
        chemistry::SmallMoleculeMiscProperties().GetMDLProperties().Has( "property"),
        false
      );

      // test that does mdl property exists does not return false negatives
      props.SetMDLProperty( "property", std::string( "1"));
      BCL_ExampleIndirectAssert
      (
        props.GetMDLProperties().Has( "property"),
        true,
        "SetMDLProperty"
      );

      // test that we get the property string back properly
      BCL_ExampleIndirectCheck( props.GetMDLProperty( "property"), "1", "GetMDLProperty");

      props.SetMDLProperty( "property", std::string( "2"));

      // check that there is really just one mdl property in the map
      chemistry::SmallMoleculeMiscProperties::const_iterator itr( props.Begin());
      BCL_ExampleIndirectCheck( props.Begin() != props.End() && ++itr == props.End(), true, "Set, Begin, End");

      // check that we can add a property
      props.SetMDLProperty( "another property", std::string( "3"));
      itr = props.Begin();
      BCL_ExampleIndirectCheck
      (
        props.Begin() != props.End() && ++itr != props.End() && ++itr == props.End(),
        true,
        "SetMDLProperty 2nd time"
      );

      // check that the new property and the original property have the correct values
      BCL_ExampleIndirectCheck
      (
        props.GetMDLProperty( "property"),
        "2",
        "SetMDLProperty(a) should not affect property b"
      );
      BCL_ExampleIndirectCheck( props.GetMDLProperty( "another property"), "3", "SetMDLProperty");

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  };

  const ExampleClass::EnumType ExampleChemistrySmallMoleculeMiscProperties::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistrySmallMoleculeMiscProperties())
  );

} // namespace bcl
