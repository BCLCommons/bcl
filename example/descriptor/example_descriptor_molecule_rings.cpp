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
#include "descriptor/bcl_descriptor_molecule_rings.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_rings.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 13, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeRings :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeRings *Clone() const
    {
      return new ExampleDescriptorMoleculeRings( *this);
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
      descriptor::MoleculeRings nrings_default;
      descriptor::MoleculeRings nrings_aromatic( chemistry::ConstitutionalBondTypeData::e_Aromatic);
      descriptor::MoleculeRings nrings_macrocycle( chemistry::ConstitutionalBondTypeData::e_Any, true);

      // copy constructor
      descriptor::MoleculeRings nrings_copy( nrings_aromatic);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( nrings_default.GetAlias(), "NRings");
      BCL_ExampleCheck( nrings_aromatic.GetAlias(), "NAromaticRings");
      BCL_ExampleCheck( nrings_macrocycle.GetAlias(), "NMacrocyclicRings");

      BCL_ExampleCheck( nrings_default.GetString(), "NRings");

    ///////////////
    // operators //
    ///////////////

      // map from filename to # of expected rings (first) and expected aromatic rings (second)
      storage::Map< std::string, storage::Vector< size_t> > filenames_to_expected_n_rings;

      // diazepam has 3 rings, 2 of which are aromatic, and no macrocycles
      filenames_to_expected_n_rings[ "diazepam.sdf"] = storage::Vector< size_t>::Create( 3, 2, 0);
      // taxol has 7 rings, 3 of which are aromatic, and no macrocycles
      filenames_to_expected_n_rings[ "taxol.sdf"] = storage::Vector< size_t>::Create( 7, 3, 0);
      // hexane has no rings
      filenames_to_expected_n_rings[ "hexane.sdf"] = storage::Vector< size_t>::Create( 0, 0, 0);
      // cyclododecane has a single non-aromatic ring which is also a macrocycle
      filenames_to_expected_n_rings[ "cyclododecane.sdf"] = storage::Vector< size_t>::Create( 1, 0, 1);

      for
      (
        storage::Map< std::string, storage::Vector< size_t> >::const_iterator
          itr_files( filenames_to_expected_n_rings.Begin()),
          itr_files_end( filenames_to_expected_n_rings.End());
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

        const chemistry::FragmentComplete &mol( ensemble.GetMolecules().FirstElement());

        BCL_ExampleIndirectCheck
        (
          nrings_default( mol).First(),
          float( itr_files->second( 0)),
          "Counting rings in " + itr_files->first
        );

        BCL_ExampleIndirectCheck
        (
          nrings_aromatic( mol).First(),
          float( itr_files->second( 1)),
          "Counting aromatic rings in " + itr_files->first
        );

        BCL_ExampleIndirectCheck
        (
          nrings_macrocycle( mol).First(),
          float( itr_files->second( 2)),
          "Counting macrocyclic rings in " + itr_files->first
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( nrings_default, descriptor::MoleculeRings()),
        true,
        "NRings I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeRings

  const ExampleClass::EnumType ExampleDescriptorMoleculeRings::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeRings())
  );

} // namespace bcl
