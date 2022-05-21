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
#include "descriptor/bcl_app_descriptor_analyze.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include "chemistry/bcl_chemistry.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_descriptor_analyze.cpp
  //!
  //! @author butkiem1
  //! @date Aug 23, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppDescriptorAnalyze :
    public ExampleInterface
  {
  public:

    ExampleAppDescriptorAnalyze *Clone() const
    {
      return new ExampleAppDescriptorAnalyze( *this);
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

      ApplicationExampleHelper descriptor_analyze_helper( app::DescriptorAnalyze::DescriptorAnalyze_Instance);

    /////////////////
    // data access //
    /////////////////

      // construct code object file for vertical and horizontal axis
      const std::string code_obj_horiz_a
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "code_horiz_a.object")
      );
      const std::string code_obj_horiz_b
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "code_horiz_b.object")
      );
      const std::string code_obj_vert_a
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "code_vert_a.object")
      );
      const std::string code_obj_vert_b
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "code_vert_b.object")
      );

      // filename for overlap matrix
      const std::string overlap_matrix_filename
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "code_obj_overlap_matrix.txt")
      );

      // filename for overlap matrix
      const std::string overlap_matrix_filename_reference
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "code_obj_overlap_matrix.reference.txt")
      );

      // reference to an arbitrary molecule file
      const std::string taxol_filename_reference
      (
        AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf")
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // create a command line
      descriptor_analyze_helper.ResetFlagsAndParameters();

      descriptor_analyze_helper.SetFlag
      (
        "code_object_files_horizontal",
        storage::Vector< std::string>::Create( code_obj_horiz_a, code_obj_horiz_b)
      );

      descriptor_analyze_helper.SetFlag
      (
        "code_object_files_vertical",
        storage::Vector< std::string>::Create( code_obj_vert_a, code_obj_vert_b)
      );
      descriptor_analyze_helper.SetFlag
      (
        "retriever",
        "SdfFile(filename=" + taxol_filename_reference + ")"
      );

      descriptor_analyze_helper.SetFlag( "output_matrix_filename", overlap_matrix_filename);

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( descriptor_analyze_helper.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( descriptor_analyze_helper.RunCommand(), 0))
      {
        // if the application ran successfully, check that files match
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( overlap_matrix_filename, overlap_matrix_filename_reference),
          true,
          "overlap matrix was generated correctly"
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end Exampleexample_name

  const ExampleClass::EnumType ExampleAppDescriptorAnalyze::s_Instance
  (
    GetExamples().AddEnum( ExampleAppDescriptorAnalyze())
  );

} // namespace bcl
