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
#include "io/bcl_io_fixed_line_width_writer.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_io_fixed_line_width_writer.cpp
  //!
  //! @author mendenjl
  //! @date Mar 1, 2013
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoFixedLineWidthWriter :
    public ExampleInterface
  {
  public:

    ExampleIoFixedLineWidthWriter *Clone() const
    {
      return new ExampleIoFixedLineWidthWriter( *this);
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

      {
        io::FixedLineWidthWriter writer( 0, 40);

        writer << "this string should be split after here to maintain 40 chars per line";

        // try word wrap
        BCL_ExampleIndirectCheck
        (
          writer.String(),
          "this string should be split after here\nto maintain 40 chars per line",
          "writer << \"this string should be split after here to maintain 40 chars per line\""
        );
      }

      {
        // try out the auto-indent continued lines extra option
        io::FixedLineWidthWriter writer_extra( 2, 40);

        writer_extra << "this string should be split after here to maintain 40 chars per line";

        // try word wrap
        BCL_ExampleIndirectCheck
        (
          writer_extra.String(),
          "this string should be split after here\n  to maintain 40 chars per line",
          "writer_extra << \"this string should be split after here to maintain 40 chars per line\""
        );
      }

      // try the set indent function
      {
        io::FixedLineWidthWriter writer( 0, 40);
        writer.SetIndent( 2);
        writer << "this string should be split after here x to maintain 40 chars per line";

        // try word wrap
        BCL_ExampleIndirectCheck
        (
          writer.String(),
          "  this string should be split after here\n  x to maintain 40 chars per line",
          "SetIndent"
        );
      }

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoFixedLineWidthWriter

  const ExampleClass::EnumType ExampleIoFixedLineWidthWriter::s_Instance
  (
    GetExamples().AddEnum( ExampleIoFixedLineWidthWriter())
  );

} // namespace bcl
