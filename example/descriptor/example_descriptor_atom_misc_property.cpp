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
#include "descriptor/bcl_descriptor_atom_misc_property.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_misc_property.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 19, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomMiscProperty :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomMiscProperty *Clone() const
    {
      return new ExampleDescriptorAtomMiscProperty( *this);
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
      // Initialize molecule
      chemistry::FragmentComplete molecule;
      // set test misc property
      molecule.StoreProperty( std::string( "example_property"), std::string( "1,3,5,7"));
      // set test misc property
      molecule.StoreProperty( std::string( "not_exist"), std::string( "3,6,9"));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // check default constructor
      descriptor::AtomMiscProperty atom_property_default;

      // check constructor with parameters
      descriptor::AtomMiscProperty atom_property( "example_property", size_t( 4));

    /////////////////
    // data access //
    /////////////////

      // check GetAlias
      BCL_ExampleCheck( atom_property.GetAlias(), "Atom_MiscProperty");

      // check GetNumberOfReturnValuesPerAtom
      BCL_ExampleCheck( atom_property.GetNormalSizeOfFeatures(), size_t( 4));

      // check Recalculate
      BCL_ExampleCheck( atom_property.CollectValuesOnEachElementOfObject( molecule).GetSize(), size_t( 0));

      // create instance of atom property
      descriptor::AtomMiscProperty atom_property_exists( "exists", size_t( 3));

      // check GetNumberOfReturnValuesPerAtom
      BCL_ExampleCheck( atom_property_exists.GetNormalSizeOfFeatures(), size_t( 3));

      // check Recalculate
      BCL_ExampleCheck( atom_property_exists.CollectValuesOnEachElementOfObject( molecule).GetSize(), size_t( 0));

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomMiscProperty

  const ExampleClass::EnumType ExampleDescriptorAtomMiscProperty::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomMiscProperty())
  );

} // namespace bcl
