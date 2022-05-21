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
#include "internal/biology/bcl_app_scatterplot_from_tables.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "example_check_macros.h"
#include "biol/bcl_biol.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_scatterplot_from_tables.cpp
  //! @brief this example tests the implementation of the app ScatterplotFromTables, which generates scatter plots from
  //! scoring tables
  //!
  //! @author fischea
  //! @date Jun 22, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppScatterplotFromTables :
     public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleAppScatterplotFromTables
    ExampleAppScatterplotFromTables *Clone() const
    {
      return new ExampleAppScatterplotFromTables( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      const app::ApplicationType app_enum_scatter( "ScatterplotFromTables");
      BCL_ExampleAssert( app_enum_scatter.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // creation of a scatter plot based on 4 models
      {
        ApplicationExampleHelper scatter_helper( app_enum_scatter);

        // get paths and filenames
        const std::string fold_table( AddExampleInputPathToFilename( e_Biology, "scatterplots/1.tbl"));
        const std::string storage_path( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), ""));
        const std::string prefix( "1_test_");
        const std::string output_gnuplot
        (
          AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), prefix + "scatter_sum_rmsd.gnuplot")
        );
        const std::string ss_path( AddExampleInputPathToFilename( e_Biology, ""));

        // test using native pool, default scores, and single stage
        scatter_helper.SetFlag( "gnuplot_input_table_filenames", fold_table);
        scatter_helper.SetFlag( "gnuplot_table_columns_x", "RMSD");
        scatter_helper.SetFlag( "gnuplot_table_columns_y", "sum");
        scatter_helper.SetFlag( "gnuplot_output_filename", output_gnuplot);
        scatter_helper.SetFlag( "gnuplot_pixels_x", "600");
        scatter_helper.SetFlag( "gnuplot_pixels_y", "400");
        scatter_helper.SetFlag( "gnuplot_set_key");
        scatter_helper.SetFlag( "gnuplot_x_label", "RMSD");
        scatter_helper.SetFlag( "gnuplot_y_label", "sum");
        scatter_helper.SetFlag( "gnuplot_min_x", "0");
        scatter_helper.SetFlag( "gnuplot_max_x", "20");
        scatter_helper.SetFlag( "gnuplot_min_y", "-60000");
        scatter_helper.SetFlag( "gnuplot_max_y", "-40000");
        scatter_helper.SetFlag( "gnuplot_tics_x", "1");
        scatter_helper.SetFlag( "gnuplot_tics_y", "2");
        scatter_helper.SetFlag( "gnuplot_series_names", "RMSD");
        scatter_helper.SetFlag( "gnuplot_title", "scatter plot test");

        // check the command line and run
        BCL_ExampleAssert( scatter_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( scatter_helper.RunCommand(), 0);

        // check that the output is unchanged

        // check the gnuplot file
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance( output_gnuplot, output_gnuplot + ".correct", 0.1),
            true,
            "Scatterplot gnuplot file was generated correctly"
          )
        )
        {
          // example check success, remove output file
          remove( output_gnuplot.c_str());
        }
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleAppScatterplotFromTables

  //! single instance of this class
  const ExampleClass::EnumType ExampleAppScatterplotFromTables::s_Instance
  (
     GetExamples().AddEnum( ExampleAppScatterplotFromTables())
  );

} // namespace bcl
