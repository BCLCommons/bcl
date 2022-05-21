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
#include "assemble/bcl_assemble_protein_with_cache_storage_file.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_factory.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_protein_with_cache_storage_file.cpp
  //!
  //! @author mendenjl
  //! @date Jan 09, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleProteinWithCacheStorageFile :
    public ExampleInterface
  {
  public:

    ExampleAssembleProteinWithCacheStorageFile *Clone() const
    {
      return new ExampleAssembleProteinWithCacheStorageFile( *this);
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
      // directory that contains a bunch of pdbs and fastas
      const std::string pdb_fasta_dir( AddExampleInputPathToFilename( e_Biology, ""));

      // create an object to read proteins, as either pdbs or fastas, from the directory
      util::Implementation
      <
        assemble::RetrieveProteinModelWithCache
      > retriever_all( "ProteinDirectory( " + pdb_fasta_dir + ")");

      // ensure that there are at least 10 proteins in the input directory
      BCL_ExampleCheck( retriever_all->GetAllKeys().GetSize() > 10, true);

      // Try to retrieve ubiquiton
      assemble::ProteinModelWithCache ubiquiton
      (
        pdb::Factory( biol::GetAAClasses().e_AA).ProteinModelFromPDBFilename
        (
          pdb_fasta_dir + "1ubi.pdb",
          storage::Map< biol::SSType, size_t>(),
          true
        ),
        false
      );
      // ensure that the sequences are the same
      BCL_ExampleCheck
      (
        retriever_all->Retrieve( "1ubi")->GetSequences()( 0)->Sequence(),
        ubiquiton.GetSequences()( 0)->Sequence()
      );

      // check sizes on other proteins in the directory
      BCL_ExampleCheck( retriever_all->Retrieve( "1ubi")->GetSize(), 76);

      // 0 because we requested Proteins, not sequences, and only the fasta file exists
      BCL_ExampleCheck( retriever_all->Retrieve( "1KSR").IsDefined(), false);

      util::SiPtr< const assemble::ProteinWithCacheStorageFile> si_ptr_storage( retriever_all.operator->());
      BCL_MessageVrb
      (
        "Here are all the pdb/fasta files found by this class: " + util::Format()( si_ptr_storage->GetAllLocations())
      );

      // create an object to read fastas from sub-directories
      util::Implementation< assemble::RetrieveProteinModelWithCache> retriever_selected_fastas_sub_directories
      (
        "SequenceDirectory( "
        + GetExamples().GetExamplePath() + s_ExampleInputFolderName + PATH_SEPARATOR + ", "
        + "key file = " + AddExampleInputPathToFilename( e_Biology, "protein_cache_storage.ls") + ","
        + "recursive=True)"
      );

      // check sizes
      BCL_ExampleAssert( retriever_selected_fastas_sub_directories.IsDefined(), true);
      BCL_ExampleCheck( retriever_selected_fastas_sub_directories->GetSize(), 4);
      BCL_ExampleCheck( retriever_selected_fastas_sub_directories->Retrieve( "1ubi")->GetSize(), 76);
      BCL_ExampleCheck( retriever_selected_fastas_sub_directories->Retrieve( "1ubi")->GetSize(), 76);
      BCL_ExampleCheck( retriever_selected_fastas_sub_directories->Retrieve( "1KSR")->GetSize(), 92);
      BCL_ExampleCheck( retriever_selected_fastas_sub_directories->Retrieve( "1B43A")->GetSize(), 340);
      BCL_ExampleCheck( retriever_selected_fastas_sub_directories->Retrieve( "cluster_0000_final")->GetSize(), 76);

      // retrieve an ensemble
      util::ShPtrVector< assemble::ProteinModelWithCache> ensemble
      (
        retriever_selected_fastas_sub_directories->RetrieveEnsemble( math::Range< size_t>( 1, 2))
      );
      BCL_ExampleIndirectCheck( ensemble.GetSize(), 2, "retrieve ensemble");
      BCL_ExampleIndirectCheck
      (
        ensemble( 0)->GetIdentification(),
        retriever_selected_fastas_sub_directories->Retrieve( "1KSR")->GetIdentification(),
        "retrieve ensemble"
      );
      BCL_ExampleIndirectCheck
      (
        ensemble( 1)->GetIdentification(),
        retriever_selected_fastas_sub_directories->Retrieve( "1B43A")->GetIdentification(),
        "retrieve ensemble"
      );

      BCL_MessageStd
      (
        "Example directory searcher " + retriever_selected_fastas_sub_directories.GetString()
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleProteinWithCacheStorageFile

  const ExampleClass::EnumType ExampleAssembleProteinWithCacheStorageFile::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleProteinWithCacheStorageFile())
  );

} // namespace bcl
