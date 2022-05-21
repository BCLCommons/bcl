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
#include "chemistry/bcl_chemistry_substituent_conformational.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_molecule_complete.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_substituent_conformational.cpp
  //!
  //! @author mendenjl
  //! @date Jan 26, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistrySubstituentConformational :
    public ExampleInterface
  {
  public:

    ExampleChemistrySubstituentConformational *Clone() const
    {
      return new ExampleChemistrySubstituentConformational( *this);
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

      // default constructor
      chemistry::SubstituentConformational substitutent_default;

      // construct from properties
      util::ShPtr< chemistry::SubstituentConformational> sp_substituent
      (
        new chemistry::SubstituentConformational()
      );

      // copy constructor
      util::ShPtr< chemistry::SubstituentConformational> substitutent_copy( sp_substituent);

      // clone
      const util::ShPtr< chemistry::SubstituentConformational> substitutent_clone( substitutent_default.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_MessageStd( "class name: " + substitutent_clone->GetClassIdentifier());
      BCL_ExampleAssert( substitutent_clone.IsDefined(), true);
      BCL_ExampleCheck( chemistry::SubstituentConformational().GetClassIdentifier(), GetStaticClassName< chemistry::SubstituentConformational>());

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // setup input stream
      io::IFStream input_sdf;

      // read in molecule
      chemistry::MoleculeComplete small_mol;
      const std::string filename_in( AddExampleInputPathToFilename( e_Chemistry, "corina_diazepam.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, filename_in);
      small_mol = sdf::Factory::MakeMolecule( input_sdf);

      storage::Vector< chemistry::SubstituentConformational> all_substituents;
      all_substituents.AllocateMemory( small_mol.GetNumberAtoms());
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms( small_mol.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        chemistry::SubstituentConformational current_substituent( *itr_atoms);
        all_substituents.PushBack( current_substituent);
      }

      //Sort SubstituentConformational vector and check that they are properly sorted in descending priority
      all_substituents.Sort( std::less< chemistry::SubstituentConformational>());

      BCL_ExampleCheck
      (
        all_substituents( 1).GetRootAtom()->GetElementType(),
        chemistry::GetElementTypes().e_Oxygen
      );

      BCL_ExampleCheck
      (
        all_substituents.FirstElement().GetRootAtom()->GetElementType(),
        chemistry::GetElementTypes().e_Chlorine
      );

      BCL_ExampleCheck
      (
        all_substituents.LastElement().GetRootAtom()->GetElementType(),
        chemistry::GetElementTypes().e_Hydrogen
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistrySubstituentConformational

  const ExampleClass::EnumType ExampleChemistrySubstituentConformational::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistrySubstituentConformational())
  );

} // namespace bcl
