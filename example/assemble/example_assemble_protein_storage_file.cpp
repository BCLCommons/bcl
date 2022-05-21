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
#include "assemble/bcl_assemble_protein_storage_file.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_protein_storage_file.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleProteinStorageFile :
    public ExampleInterface
  {
  public:

    ExampleAssembleProteinStorageFile *Clone() const
    { return new ExampleAssembleProteinStorageFile( *this);}

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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      //read pdbfile
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AACaCb, ssetype_min_size)
      );

      // initializer
      io::Directory initializer( AddExampleOutputPathToFilename( protein_model, "run_aa"));

      // remove the storage recursive in case it was not properly deleted last run
      initializer.Remove( true);

      // sources
      const std::string source_a( "source1");
      const std::string source_b( "source2");

    //////////////////
    // construction //
    //////////////////

      // construct
      assemble::ProteinStorageFile protein_storage_attach( initializer.GetPath(), assemble::ProteinStorageFile::e_Attach);
      assemble::ProteinStorageFile protein_storage_overwrite( initializer.GetPath(), assemble::ProteinStorageFile::e_Overwrite);

      // clone
      util::ShPtr< assemble::ProteinStorageFile> sp_storage( protein_storage_attach.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( sp_storage->GetClassIdentifier(), GetStaticClassName( protein_storage_attach));

      // Test that attach and overwrite both can manage the same storage

      // get initializer
      BCL_ExampleCheck( protein_storage_overwrite.GetInitializer(), initializer.GetPath());
      BCL_ExampleCheck( protein_storage_attach.GetInitializer(), initializer.GetPath());

      // get all keys
      BCL_ExampleCheck( protein_storage_attach.GetAllKeys( source_a).GetSize(), 0);
      BCL_ExampleCheck( protein_storage_attach.GetSize( source_a), 0);

    ////////////////
    // operations //
    ////////////////

      // storage
      const std::string store_a1_key( protein_storage_overwrite.Store( protein_model, source_a));
      const std::string store_a2_key( protein_storage_attach.Store( protein_model, source_a));

      BCL_ExampleIndirectCheck( store_a1_key, "_0000", "Store attach");
      BCL_ExampleIndirectCheck( store_a2_key, "_0001", "Store create");

      // store using a key
      BCL_ExampleCheck( protein_storage_attach.Store( protein_model, source_a, store_a1_key), true);
      const std::string store_a3_key( "_0002");
      BCL_ExampleCheck( protein_storage_attach.Store( protein_model, source_a, store_a3_key), true);

      // store ensemble
      util::ShPtr< assemble::ProteinModel> sp_model( protein_model.Clone());
      util::ShPtrList< assemble::ProteinModel> ensemble( 2, sp_model);
      storage::Vector< std::string> store_a4a5_keys( protein_storage_attach.Store( ensemble, source_a));
      BCL_ExampleIndirectCheck( store_a4a5_keys.GetSize(), 2, "store ensemble");
      BCL_ExampleIndirectCheck( store_a4a5_keys( 0), "_0004", "store ensemble");
      BCL_ExampleIndirectCheck( store_a4a5_keys( 1), "_0005", "store ensemble");

      // check that sourceb is still empty
      BCL_ExampleIndirectCheck( protein_storage_attach.GetAllKeys( source_b).GetSize(), 0, "source b empty");

      // retrieve
      util::ShPtr< assemble::ProteinModel> sp_model1( protein_storage_attach.Retrieve( source_a, store_a1_key));
      util::ShPtr< assemble::ProteinModel> sp_model2( protein_storage_attach.Retrieve( source_b, "000007"));
      BCL_ExampleIndirectCheck( sp_model1.IsDefined(), true, "retrieve by source and key1");
      BCL_ExampleIndirectCheck( sp_model2.IsDefined(), false, "retrieve by source and wrong key");

      // retrieve ensemble
      BCL_ExampleIndirectCheck( protein_storage_attach.RetrieveEnsemble( source_a).GetSize(), 6, "retrieve ensemble for source");
      BCL_ExampleIndirectCheck( protein_storage_attach.RetrieveEnsemble( source_b).GetSize(), 0, "retrieve ensemble for source");

      util::ShPtrList< assemble::ProteinModel> ensemble_retrieved( protein_storage_attach.RetrieveEnsemble( source_a, store_a4a5_keys));
      BCL_ExampleIndirectCheck( ensemble_retrieved.GetSize(), store_a4a5_keys.GetSize(), "retrieve ensemble for source from keys");
      BCL_ExampleIndirectCheck( protein_storage_attach.RetrieveEnsemble( source_b, store_a4a5_keys).IsEmpty(), true, "retrieve by keys from wrong source");

      // retrieve by range
      BCL_ExampleIndirectCheck( protein_storage_attach.RetrieveEnsemble( source_a, math::Range< size_t>( 0, 1)).GetSize(), 2, "retrieve by range");

    //////////////////////
    // input and output //
    //////////////////////

      // read write
      WriteBCLObject( protein_storage_attach);
      assemble::ProteinStorageFile protein_storage_read;
      ReadBCLObject( protein_storage_read);

      BCL_ExampleIndirectCheck( protein_storage_attach.GetInitializer(), protein_storage_read.GetInitializer(), "read write");
      BCL_ExampleIndirectCheck( protein_storage_attach.GetAllKeys( source_a).GetSize(), protein_storage_read.GetAllKeys( source_a).GetSize(), "read write");

    /////////////
    // cleanup //
    /////////////

      // remove the storage recursive
      initializer.Remove( true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // ExampleAssembleProteinStorageFile

  const ExampleClass::EnumType ExampleAssembleProteinStorageFile::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleProteinStorageFile())
  );

} // namespace bcl
