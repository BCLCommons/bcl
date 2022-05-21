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
#include "biol/bcl_biol_atom_types.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_atom_types.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAtomTypes :
    public ExampleInterface
  {
  public:

    ExampleBiolAtomTypes *Clone() const
    { return new ExampleBiolAtomTypes( *this);}

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
      BCL_MessageStd( "name of the first atom type enum = " + biol::AtomType( size_t( 0)).GetName());
      BCL_MessageStd( "atom name of the second  atom type = " + util::Format()( biol::AtomType( 1)->GetName()));
      BCL_MessageStd( "Element Types as enum of the third atom type = "      + util::Format()( biol::AtomType(  3)->GetElementType()));
      BCL_MessageStd( "Atom Type as enum of the fourth atom type = " + util::Format()( biol::AtomType( 3)));
      BCL_MessageStd( "Atom Type as enum from atom name = " + util::Format()( biol::AtomType( "CB")));

      // check that ca and cb are properly stored in that vector of atom types
      const storage::Vector< biol::AtomType> ca_cb( biol::GetAtomTypes().GetCACB().Begin(), biol::GetAtomTypes().GetCACB().End());

      BCL_Example_Check
      (
        ca_cb.GetSize() == 2,
        "number of atomtypes in ca_cb should be two but is " + util::Format()( ca_cb.GetSize())
      );
      BCL_Example_Check
      (
        ca_cb( 0) == biol::GetAtomTypes().CA && ca_cb( 1) == biol::GetAtomTypes().CB,
        "the first and second element in the ca cb atom types should be CA and CB but instead are "
        + util::Format()( ca_cb( 0)) + " and " + util::Format()( ca_cb( 1))
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAtomTypes

  const ExampleClass::EnumType ExampleBiolAtomTypes::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAtomTypes())
  );

} // namespace bcl

