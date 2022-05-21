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
#include "model/bcl_model_meta_data_storage_file.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_meta_data_storage_file.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author butkiem1
  //! @date Jan 26, 2012
  //! @remarks status incomplete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelMetaDataStorageFile :
    public ExampleInterface
  {
  public:

    ExampleModelMetaDataStorageFile *Clone() const
    {
      return new ExampleModelMetaDataStorageFile( *this);
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
      // sample data for storage object
      const float result( 0.123);
      const util::ObjectDataLabel descriptors( "Combine()");
      const util::ObjectDataLabel method_name( "Combine()");
      const util::ObjectDataLabel obj_function( "Combine()");

      const descriptor::Dataset exp_pred;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      model::MetaDataStorageFile storage_default;

      //! construct a storage interface implementation
      util::Implementation< model::MetaDataStorageInterface> storage
      (
        "File( directory=" + AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "dir_storage_meta") + ")"
      );

    /////////////////
    // data access //
    /////////////////

      // store model
      BCL_ExampleCheck
      (
        storage->Store( result, descriptors, method_name, obj_function),
        "000000"
      );

      // store model
      BCL_ExampleCheck
      (
        storage->Store( result, descriptors, method_name, obj_function, "000099"),
        "000099"
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

      // write bcl object
      WriteBCLObject( storage_default);

      //! default object
      model::MetaDataStorageFile storage_read;
      // read bcl object
      ReadBCLObject( storage_read);

      // clean up created directories
      io::Directory dir( AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "dir_storage_meta"));
      dir.Remove( true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelMetaDataStorageFile

  const ExampleClass::EnumType ExampleModelMetaDataStorageFile::s_Instance
  (
    GetExamples().AddEnum( ExampleModelMetaDataStorageFile())
  );

} // namespace bcl
