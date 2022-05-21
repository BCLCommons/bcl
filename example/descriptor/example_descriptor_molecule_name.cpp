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
#include "descriptor/bcl_descriptor_molecule_name.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_name.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 17, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeName :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeName *Clone() const
    {
      return new ExampleDescriptorMoleculeName( *this);
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

      // create a property that retrieves the name
      descriptor::CheminfoID name_retriever( "Name");

      const std::string test_name( " test name");
      chemistry::FragmentComplete fragment_with_name( chemistry::AtomVector< chemistry::AtomComplete>(), test_name);

      // ensure that each is defined
      BCL_ExampleAssert( name_retriever.IsDefined(), true);

      // check the descriptor string
      BCL_ExampleCheck( name_retriever.GetString(), "Name(size=240)");

      // try each out on a small molecule
      BCL_ExampleCheck
      (
        std::string ( name_retriever->SumOverObject( fragment_with_name).Begin()),
        util::TrimString( fragment_with_name.GetName())
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeName

  const ExampleClass::EnumType ExampleDescriptorMoleculeName::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeName())
  );

} // namespace bcl
