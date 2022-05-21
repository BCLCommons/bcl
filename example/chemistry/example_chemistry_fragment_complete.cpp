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
#include "chemistry/bcl_chemistry_fragment_complete.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_complete.cpp
  //! @details Tests ChemistryFragmentComplete class which contains small molecule configuration data
  //!
  //! @author kothiwsk
  //! @date
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentComplete :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentComplete *Clone() const
    {
      return new ExampleChemistryFragmentComplete( *this);
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

      // construct from conformation interface
      io::IFStream input_sdf;
      const std::string hexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, hexane_filename);

      // load information into small_mol_conformation
      chemistry::FragmentComplete small_mol_conformation;
      small_mol_conformation
        = sdf::FragmentFactory::MakeFragment( input_sdf);
      // close the input stream
      io::File::CloseClearFStream( input_sdf);

    /////////////////
    // data access //
    /////////////////

      // check GetAtomsIterator
      size_t atom_index = 0;
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface>
        itr( small_mol_conformation.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        if( itr->GetAtomType() != chemistry::GetAtomTypes().C_TeTeTeTe)
        {
          break;
        }
        ++atom_index;
      }
      // there should be 6 atoms of C_TeTeTeTe type
      BCL_ExampleIndirectCheck( atom_index, 6, "Conformational iterator succeeds to iterate over underlying atom vector");

      // test SetName and GetName
      small_mol_conformation.SetName( "hexane");
      BCL_ExampleCheck( small_mol_conformation.GetName(), "hexane");

      // test GetBonds
      BCL_ExampleCheck( small_mol_conformation.GetNumberBonds(), 6);

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

  }; //end ExampleSmallMolecule

  const ExampleClass::EnumType ExampleChemistryFragmentComplete::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentComplete())
  );

} // namespace bcl
