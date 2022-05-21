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
#include "descriptor/bcl_descriptor_molecule_rotatable_bonds.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_rotatable_bonds.cpp
  //!
  //! @author kothiwsk, mendenjl
  //! @date Feb 13, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeRotatableBonds :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeRotatableBonds *Clone() const
    {
      return new ExampleDescriptorMoleculeRotatableBonds( *this);
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
      descriptor::MoleculeRotatableBonds nrot;

      // copy constructor
      descriptor::MoleculeRotatableBonds nrot_copy( nrot);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( nrot.GetAlias(), "NRotBond");

      BCL_ExampleCheck( nrot.GetString(), "NRotBond");

    ///////////////
    // operators //
    ///////////////

      storage::Map< std::string, size_t> filenames_to_expected_n_rotatable_bonds;
      // the values used for expected #s of rotatable bonds can also be found in sdf files downloaded from pubchem,
      // under the descriptor called PUBCHEM_CACTVS_ROTATABLE_BONDS

      // diazepam has 1 rotatable bond (the bond between the benzene and fused ring system)
      filenames_to_expected_n_rotatable_bonds[ "diazepam.sdf"] = 1;
      filenames_to_expected_n_rotatable_bonds[ "taxol.sdf"] = 14; // taxol has 14 rotatable bonds (1 amide bond ignored)
      filenames_to_expected_n_rotatable_bonds[ "hexane.sdf"] = 3; // the 3 internal CC bonds of hexane can be rotated

      for
      (
        storage::Map< std::string, size_t>::const_iterator
          itr_files( filenames_to_expected_n_rotatable_bonds.Begin()),
          itr_files_end( filenames_to_expected_n_rotatable_bonds.End());
        itr_files != itr_files_end;
        ++itr_files
      )
      {
        // create input stream for reading a smallmolecule ensemble
        io::IFStream input;
        BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, itr_files->first));
        // read in ensemble
        chemistry::FragmentEnsemble ensemble( input);
        // close stream
        io::File::CloseClearFStream( input);

        const chemistry::FragmentComplete &mol( *ensemble.Begin());

        BCL_ExampleIndirectCheck
        (
          nrot( mol).First(),
          float( itr_files->second),
          "Counting rotatable bonds in " + itr_files->first
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( nrot, nrot_copy),
        true,
        "NRotBond I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeRotatableBonds

  const ExampleClass::EnumType ExampleDescriptorMoleculeRotatableBonds::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeRotatableBonds())
  );

} // namespace bcl
