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
#include "model/bcl_model_interface_retrieve_from_file.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_file.h"
#include "model/bcl_model_decision_tree.h"
#include "model/bcl_model_store_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_interface_retrieve_from_file.cpp
  //! @brief check retrieval methods of model::InterfaceRetrieveFromFile using the file system.
  //!
  //! @author butkiem1
  //! @date Feb 24, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelInterfaceRetrieveFromFile :
    public ExampleInterface
  {
  public:

    ExampleModelInterfaceRetrieveFromFile *Clone() const
    {
      return new ExampleModelInterfaceRetrieveFromFile( *this);
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
      util::ShPtr< model::Interface> model( new model::DecisionTree());
      const float result( 0.0);

      const util::ObjectDataLabel descriptors( "Combine()");
      const util::ObjectDataLabel method_name( "Combine()");
      const util::ObjectDataLabel obj_function( "Combine()");
      const util::ObjectDataLabel result_label( "PredictedActivity");

      //! construct a storage interface implementation
      util::Implementation< model::StoreInterface> storage_write
      (
        "File( directory=" + AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "dir_retrieve_model") + ")"
      );

      storage_write->StoreResultDescriptor( result_label);

      // store model
      storage_write->Store( model, result, descriptors, method_name, obj_function, model::CrossValidationInfo());

      // store model
      storage_write->Store( model, result, descriptors, method_name, obj_function, model::CrossValidationInfo(), "000099");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      model::InterfaceRetrieveFromFile storage_default;

      //! construct a retrieve interface implementation
      util::Implementation< model::RetrieveInterface> storage
      (
        "File( directory=" + AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "dir_retrieve_model") + ")"
      );

    /////////////////
    // data access //
    /////////////////

      // retrieve the results descriptor
      BCL_ExampleCheck( storage->RetrieveResultDescriptor(), result_label);

      // retrieve all models
      BCL_ExampleCheck( storage->RetrieveEnsemble().GetSize(), size_t( 2));

      // retrieve by id
      BCL_ExampleCheck( storage->RetrieveEnsemble( storage::Vector< std::string>::Create( "98", "99", "100")).GetSize(), size_t( 1));

      // retrieve by id
      BCL_ExampleCheck( storage->RetrieveEnsemble( storage::Vector< std::string>::Create( "0", "99")).GetSize(), size_t( 2));

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // setup storage default identically to the storage implementation
      BCL_ExampleAssert( storage_default.TryRead( storage.GetLabel(), util::GetLogger()), true);
      WriteBCLObject( storage_default);

      //! default object
      model::InterfaceRetrieveFromFile storage_read;
      // read bcl object
      ReadBCLObject( storage_read);

      // clean up created directories
      io::Directory dir
      (
        AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "dir_retrieve_model")
      );
      dir.Remove( true);

      util::Implementation< model::RetrieveInterface> mergetester
      (
        "File( directory=" + AddExampleInputPathToFilename( e_Model, "testmerge") + ",prefix=model)"
      );
      util::ShPtr< descriptor::Dataset> sp_merged( mergetester->ReadMergedIndependentPredictions());
      const std::string merge_file( AddExampleInputPathToFilename( e_Model, "testmerge/model.independent"));
      model::CrossValidationInfo::WritePredictions
      (
        sp_merged->GetFeatures(),
        sp_merged->GetResults(),
        sp_merged->GetIds(),
        merge_file
      );
      if
      (
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance
          (
            merge_file,
            merge_file + ".correct",
            0.00001
          ),
          true,
          "ReadMergedIndependentPredictions"
        )
      )
      {
        remove( merge_file.c_str());
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelInterfaceRetrieveFromFile

  const ExampleClass::EnumType ExampleModelInterfaceRetrieveFromFile::s_Instance
  (
    GetExamples().AddEnum( ExampleModelInterfaceRetrieveFromFile())
  );

} // namespace bcl
