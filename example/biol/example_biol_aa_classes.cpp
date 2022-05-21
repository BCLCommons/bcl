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
#include "biol/bcl_biol_aa_classes.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_classes.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAAClasses :
    public ExampleInterface
  {
  public:

    ExampleBiolAAClasses *Clone() const
    { return new ExampleBiolAAClasses( *this);}

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

      // construct AAClass from each class type
      BCL_MessageStd( "Constructing AAClass objects from all known AAClasses");

      biol::AAClass aaclass_aa( biol::GetAAClasses().e_AA);
      BCL_ExampleAssert( biol::AAClass( biol::GetAAClasses().e_AA), biol::GetAAClasses().e_AA);

      biol::AAClass aaclass_aa_ca_cb( biol::GetAAClasses().e_AACaCb);
      BCL_ExampleAssert( biol::AAClass( biol::GetAAClasses().e_AACaCb), biol::GetAAClasses().e_AACaCb);

      biol::AAClass aaclass_aa_back_bone( biol::GetAAClasses().e_AABackBone);
      BCL_ExampleAssert( biol::AAClass( biol::GetAAClasses().e_AABackBone), biol::GetAAClasses().e_AABackBone);

      biol::AAClass aaclass_aa_complete( biol::GetAAClasses().e_AAComplete);
      BCL_ExampleAssert( biol::AAClass( biol::GetAAClasses().e_AAComplete), biol::GetAAClasses().e_AAComplete);

      // construct undefined AAClass
      BCL_MessageStd( "Constructing an undefined AAClass");
      biol::AAClass aaclass_undefined( util::GetUndefined< biol::AAClass>());
      BCL_ExampleCheck( biol::AAClass( util::GetUndefined< biol::AAClass>()).IsDefined(), false);

      // use copy constructor
      BCL_MessageStd( "Calling copy constructor");
      biol::AAClass aaclass_aa_back_bone_copy( aaclass_aa_back_bone);
      BCL_ExampleCheck( biol::AAClass( aaclass_aa_back_bone), aaclass_aa_back_bone);

    /////////////////
    // data access //
    /////////////////

      // example check
      BCL_ExampleCheck( biol::GetAAClasses().GetClassIdentifier(), "bcl::biol::AAClasses");

      // test total number of AAClasses
      BCL_ExampleCheck( biol::GetAAClasses().GetEnumCount(), 4);

    ////////
    // AA //
    ////////

      // display atom types for AA
      BCL_MessageStd
      (
        "Atom Types for AA\n" +
        util::Format()( ( ( *biol::GetAAClasses().e_AA)->GetTypesOfAtoms()))
      );

      // display number of atom types for AA
      BCL_MessageStd
      (
        "The number of atom types AA has : " +
        util::Format()( ( *biol::GetAAClasses().e_AA)->GetTypesOfAtoms().GetSize())
      );

      // test number of atom types for AA
      BCL_ExampleCheck( ( *biol::GetAAClasses().e_AA)->GetTypesOfAtoms().GetSize(), 0);

    ////////////
    // AACaCb //
    ////////////

      // display atom types for AACaCb
      BCL_MessageStd
      (
        "Atom Types for AACaCb\n" + util::Format()( ( *biol::GetAAClasses().e_AACaCb)->GetTypesOfAtoms())
      );

      // test number of atom types for AACaCb
      BCL_ExampleCheck( ( *biol::GetAAClasses().e_AACaCb)->GetTypesOfAtoms().GetSize(), 2);

    ////////////////
    // AABackBone //
    ////////////////

      // display atom types for AABackBone
      BCL_MessageStd
      (
        "Atom Types for AABackBone\n" +
        util::Format()( ( *biol::GetAAClasses().e_AABackBone)->GetTypesOfAtoms())
      );

      // test number of atom types for AABackBone
      BCL_ExampleCheck( ( *biol::GetAAClasses().e_AABackBone)->GetTypesOfAtoms().GetSize(), 6);

    ////////////////
    // AAComplete //
    ////////////////

      // display atom types for AAComplete
      BCL_MessageStd
      (
        "Atom Types for AAComplete\n" +
        util::Format()( ( *biol::GetAAClasses().e_AAComplete)->GetTypesOfAtoms())
      );

      // test number of atom types for AAComplete
      BCL_ExampleCheck( ( *biol::GetAAClasses().e_AAComplete)->GetTypesOfAtoms().GetSize(), 0);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAAClasses

  const ExampleClass::EnumType ExampleBiolAAClasses::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAAClasses())
  );

} // namespace bcl

