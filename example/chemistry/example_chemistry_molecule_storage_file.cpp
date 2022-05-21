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
#include "chemistry/bcl_chemistry_molecule_storage_file.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_molecule_ensemble.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_factory.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_molecule_storage_file.cpp
  //!
  //! @author mendenjl
  //! @date Mar 09, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryMoleculeStorageFile :
    public ExampleInterface
  {
  public:

    ExampleChemistryMoleculeStorageFile *Clone() const
    {
      return new ExampleChemistryMoleculeStorageFile( *this);
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
      // example sdf filename
      const std::string diazepam_filename( AddExampleInputPathToFilename( e_Chemistry, "diazepam.sdf"));
      const std::string taxol_filename( AddExampleInputPathToFilename( e_Chemistry, "taxol_no_h.sdf"));

      //read sdf file

      // initialize input stream
      io::IFStream read;
      // input stream for diazepam
      BCL_ExampleMustOpenInputFile( read, diazepam_filename);
      // create shptr on ensemble and save molecules from sdf file
      chemistry::MoleculeComplete diazepam( sdf::Factory::MakeMolecule( read));
      io::File::CloseClearFStream( read);
      // next, read in taxol
      BCL_ExampleMustOpenInputFile( read, taxol_filename);
      chemistry::MoleculeComplete taxol( sdf::Factory::MakeMolecule( read));
      io::File::CloseClearFStream( read);

      // initializer
      std::string initializer( AddExampleOutputPathToFilename( diazepam, "smallmol_storage_file.sdf"));

    //////////////////
    // construction //
    //////////////////

      // create implementations
      util::Implementation< io::StoreInterface< chemistry::ConformationInterface> >
        molecule_storage( "File(" + initializer + ")");
      util::Implementation< io::RetrieveInterface< chemistry::MoleculeComplete, chemistry::MoleculeEnsemble> >
        molecule_retrieve( "File(" + initializer + ")");

    /////////////////
    // data access //
    /////////////////

      // check that the storage is initially empty
      BCL_ExampleCheck( molecule_retrieve->GetAllKeys().GetSize(), 0);
      BCL_ExampleCheck( molecule_retrieve->GetSize(), 0);

    ////////////////
    // operations //
    ////////////////

      // storage
      BCL_MessageStd
      (
        "Before calling store the first time: " + util::Format()( diazepam.GetBondInfo())
      );
      const std::string store_a1_key( molecule_storage->Store( diazepam));
      const std::string store_a2_key( molecule_storage->Store( taxol));

      BCL_MessageStd
      (
        "After calling store the first time: " + util::Format()( diazepam.GetBondInfo())
      );
      BCL_ExampleIndirectCheck( store_a1_key, "0", "Stored key a1");
      BCL_ExampleIndirectCheck
      (
        molecule_storage->Store( diazepam),
        store_a1_key,
        "Storing a molecule already in the storage returns the key of that molecule"
      );

      // check get key size
      util::Implementation< io::RetrieveInterface< chemistry::MoleculeComplete, chemistry::MoleculeEnsemble> >
        molecule_retrieve_new( "File(" + initializer + ")");
      BCL_ExampleCheck( molecule_retrieve_new->GetKeySize( store_a1_key), diazepam.GetNumberAtoms());
      BCL_ExampleCheck( molecule_retrieve_new->GetKeySize( store_a2_key), taxol.GetNumberAtoms());

      BCL_ExampleIndirectCheck( store_a2_key, "1", "storing a new molecule should increment the key");

      // store ensemble; translate the molecules some to change their positions so that they are added to the storage
      // again
      chemistry::MoleculeComplete diazepam_translated( diazepam);
      diazepam_translated.Translate( linal::Vector3D().SetRandomTranslation( 1.0));
      chemistry::MoleculeComplete taxol_translated( taxol);
      taxol_translated.Translate( linal::Vector3D().SetRandomTranslation( 1.0));

      storage::Vector< std::string> store_a4a5_keys;
      store_a4a5_keys.PushBack( molecule_storage->Store( diazepam_translated));
      store_a4a5_keys.PushBack( molecule_storage->Store( taxol_translated));
      BCL_ExampleIndirectCheck( store_a4a5_keys.GetSize(), 2, "store ensemble");
      BCL_ExampleIndirectCheck( store_a4a5_keys( 0), "2", "store ensemble");
      BCL_ExampleIndirectCheck( store_a4a5_keys( 1), "3", "store ensemble");

      // check that sourceb is still empty
      BCL_ExampleIndirectCheck( molecule_retrieve->GetAllKeys().GetSize(), 4, "get all keys");

      // retrieve
      chemistry::MoleculeComplete sp_model1( molecule_retrieve->Retrieve( store_a1_key));
      BCL_ExampleIndirectCheck( sp_model1.GetNumberAtoms(), 33, "retrieve by source and key1");

      chemistry::MoleculeEnsemble ensemble_retrieved( molecule_retrieve->RetrieveEnsemble( store_a4a5_keys));
      BCL_ExampleIndirectCheck( ensemble_retrieved.GetSize(), store_a4a5_keys.GetSize(), "retrieve ensemble for source from keys");

      // retrieve by range
      BCL_ExampleIndirectCheck
      (
        molecule_retrieve->RetrieveEnsemble( math::Range< size_t>( 0, 1)).GetSize(),
        2,
        "retrieve by range"
      );

    //////////////////////
    // input and output //
    //////////////////////

    /////////////
    // cleanup //
    /////////////

      // clean up generated storage file
      io::DirectoryEntry( initializer).Remove();

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryMoleculeStorageFile

  const ExampleClass::EnumType ExampleChemistryMoleculeStorageFile::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryMoleculeStorageFile())
  );

} // namespace bcl
