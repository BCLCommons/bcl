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
#include "biol/bcl_biol_atom.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_atom.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAtom :
    public ExampleInterface
  {
  public:

    ExampleBiolAtom *Clone() const
    { return new ExampleBiolAtom( *this);}

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

      const linal::Vector3D coordinates_a( 3.0, 4.0, 5.0);
      const linal::Vector3D coordinates_b( -3.0, -4.0, -5.0);

      biol::Atom atom_a( biol::GetAtomTypes().CA, 10, 0.532);
      biol::Atom atom_b( biol::AtomType( "CA"), 11, 0.532);
      biol::Atom atom_c( coordinates_a, biol::AtomType( "N"), 7, 0.242);
      biol::Atom atom_d( coordinates_b, biol::AtomType( "CB"), 8, 0.532);

      BCL_MessageStd( "name of the atom A= "  + atom_a.GetType().GetName());
      BCL_MessageStd( "name of the atom B= "  + atom_b.GetType().GetName());
      BCL_MessageStd( "pdb id of the atom A= "  + util::Format()( atom_a.GetPdbID()));
      BCL_MessageStd( "pdb id of the atom B= "  + util::Format()( atom_b.GetPdbID()));
      BCL_MessageStd( "b factor id of the atom A= "  + util::Format()( atom_a.GetBFactor()));
      BCL_MessageStd( "b factor id of the atom B= "  + util::Format()( atom_b.GetBFactor()));

      BCL_MessageStd( "name of the atom C= "  + atom_c.GetType().GetName());
      BCL_MessageStd( "name of the atom D= "  + atom_d.GetType().GetName());
      BCL_MessageStd( "coordinates of the atom C= "  + util::Format()( atom_c.GetCoordinates()));
      BCL_MessageStd( "coordinates of the atom D= "  + util::Format()( atom_d.GetCoordinates()));
      BCL_MessageStd( "element type of the atom C= " + util::Format()( atom_c.GetType()->GetElementType()));
      BCL_MessageStd( "element type of the atom D= " + util::Format()( atom_d.GetType()->GetElementType()));

      BCL_MessageStd
      (
        "distance between atoms C and D= " + util::Format()( biol::Distance( atom_c, atom_d))
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAA

  const ExampleClass::EnumType ExampleBiolAtom::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAtom())
  );

} // namespace bcl
