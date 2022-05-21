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
#include "model/bcl_model_interface_store_in_file.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "model/bcl_model_decision_tree.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_interface_store_in_file.cpp
  //! @brief check public storage and retrieval methods of model::InterfaceStoreInFile using the file system.
  //!
  //! @author butkiem1
  //! @date Feb 25, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelInterfaceStoreInFile :
    public ExampleInterface
  {
  public:

    ExampleModelInterfaceStoreInFile *Clone() const
    {
      return new ExampleModelInterfaceStoreInFile( *this);
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

      const storage::Pair< util::ShPtr< model::Interface>, float> model_result
      (
        util::ShPtr< model::Interface>( new model::DecisionTree()),
        float( 0.0)
      );

      const util::ShPtr< model::Interface> &model( model_result.First());
      const float &result( model_result.Second());

      const util::ObjectDataLabel descriptors( "Combine()");
      const util::ObjectDataLabel method_name( "Combine()");
      const util::ObjectDataLabel obj_function( "Combine()");

      const descriptor::Dataset exp_pred;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      model::InterfaceStoreInFile storage_default;

      //! construct a storage interface implementation
      util::Implementation< model::StoreInterface> storage
      (
        "File( directory=" + AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "dir_storage_model") + ")"
      );

    /////////////////
    // data access //
    /////////////////

      model::CrossValidationInfo cv_info;
      // store model
      BCL_ExampleCheck
      (
        storage->Store( model, result, descriptors, method_name, obj_function, cv_info),
        "000000"
      );

      // store model
      BCL_ExampleCheck
      (
        storage->Store( model, result, descriptors, method_name, obj_function, cv_info, "000099"),
        "000099"
      );

      //! construct a storage interface implementation
      storage = util::Implementation< model::StoreInterface>
      (
        "File( directory=" + AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "dir_storage_model") + ""
        ", key=777)"
      );

      // store model
      BCL_ExampleCheck
      (
        storage->Store( model, result, descriptors, method_name, obj_function, cv_info),
        "000777"
      );

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
      model::InterfaceStoreInFile storage_read;
      // read bcl object
      ReadBCLObject( storage_read);

      // clean up created directories
      io::Directory dir
      (
        AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "dir_storage_model")
      );
      dir.Remove( true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelInterfaceStoreInFile

  const ExampleClass::EnumType ExampleModelInterfaceStoreInFile::s_Instance
  (
    GetExamples().AddEnum( ExampleModelInterfaceStoreInFile())
  );

} // namespace bcl
