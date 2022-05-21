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
#include "util/bcl_util_object_data_label.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_object_data_label.cpp
  //!
  //! @author mendenjl
  //! @date August 01, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilObjectDataLabel :
    public ExampleInterface
  {
  public:

    ExampleUtilObjectDataLabel *Clone() const
    {
      return new ExampleUtilObjectDataLabel( *this);
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

      // construct a several valid property data labels
      util::ObjectDataLabel empty( "");
      util::ObjectDataLabel value( "0.1");
      util::ObjectDataLabel call_square( "square(            0.1     )");
      util::ObjectDataLabel call_equal( "equal( apples, oranges)");
      util::ObjectDataLabel nested_equal( "nest(nest(equal( apples, oranges)))");
      util::ObjectDataLabel nameless( "(equal( apples, oranges))");
      util::ObjectDataLabel nameless2( "(2)");
      util::ObjectDataLabel name_with_delimiter( "\"C=O\"=\"BondType=Aromatic\"(Single)");

      // construct some property data labels with quotes
      util::ObjectDataLabel nameless2_with_quotes( "(\"2\")");
      util::ObjectDataLabel quoted_stuff( "\"stuff(2)\"");

      // load in a label
      io::IFStream input_label;
      BCL_ExampleMustOpenInputFile( input_label, AddExampleInputPathToFilename( e_Chemistry, "smallmolecule_qsar_storage_code.object"));
      util::ObjectDataLabel label_from_stream( input_label);
      io::File::CloseClearFStream( input_label);

    /////////////////
    // data access //
    /////////////////

      // check that there was no string in empty
      BCL_ExampleCheck( empty.GetValue().size(), 0);

      // check that there was no string in nameless
      BCL_ExampleCheck( nameless.GetValue().size(), 0);

      // check that there was no string in nameless
      BCL_ExampleCheck( nameless2.GetValue().size(), 0);

      // check that there were no arguments to empty
      BCL_ExampleCheck( empty.GetNumberArguments(), 0);

      // check that there was a string to value
      BCL_ExampleCheck( value.GetValue(), "0.1");

      // check that there were no arguments tp value
      BCL_ExampleCheck( value.GetNumberArguments(), 0);

      // check that there was a string to square
      BCL_ExampleCheck( call_square.GetValue(), "square");

      // check that there was 1 argument to square
      BCL_ExampleCheck( call_square.GetNumberArguments(), 1);

      if( call_square.GetNumberArguments() == 1)
      {
        // check that the one argument was 0.1
        BCL_ExampleCheck( call_square.GetArgument( 0).GetValue(), "0.1");
      }

      // check that there were 2 arguments to call_equal
      BCL_ExampleCheck( call_equal.GetNumberArguments(), 2);

      if( call_equal.GetNumberArguments() == 2)
      {
        // check that the 1st argument was apples
        BCL_ExampleCheck( call_equal.GetArgument( 0).GetValue(), "apples");

        // check that the 1st argument was apples
        BCL_ExampleCheck( call_equal.GetArgument( 1).GetValue(), "oranges");
      }

      // check that there was 1 argument to nameless
      BCL_ExampleCheck( nameless.GetNumberArguments(), 1);

      // check that there were 1 arguments to nameless
      if( nameless.GetNumberArguments() == 1)
      {
        // and that it is written out the same as call_equal, even though they have different spacing
        BCL_ExampleCheck
        (
          nameless.GetArgument( 0).ToString(),
          call_equal.ToString()
        );
      }

      // check that there were 1 arguments to nameless
      if( BCL_ExampleCheck( nameless2.GetNumberArguments(), 1))
      {
        // and that it is written out the same as call_equal, even though they have different spacing
        BCL_ExampleCheck
        (
          nameless2.GetArgument( 0).ToString(),
          "2"
        );
      }

      // check that there were 1 arguments to nameless2_with_quotes
      if( BCL_ExampleCheck( nameless2_with_quotes.GetNumberArguments(), 1))
      {
        // and that it is written out the same as call_equal, even though they have different spacing
        BCL_ExampleCheck
        (
          nameless2_with_quotes.GetArgument( 0).ToString(),
          "2"
        );
      }

      // make sure there was one argument
      BCL_ExampleCheck( label_from_stream.GetNumberArguments(), 1);
      // check the reassembled label
      BCL_ExampleCheck
      (
        label_from_stream.ToString(),
        "Combine(2DA(steps=11,property=Atom_Identity,normalized=1))"
      );
      // check that comparison respects sequence order but not map/object parameter order
      BCL_ExampleCheck
      (
        util::ObjectDataLabel( "sequence(2DA(property=Atom_Identity,normalized=1,steps=11),1DA(alpha,omega))"),
        util::ObjectDataLabel( "sequence(2DA(steps=11,property=Atom_Identity,normalized=1),1DA(alpha,omega))")
      );
      BCL_ExampleIndirectCheck
      (
        util::ObjectDataLabel( "1DA(alpha,omega)") ==
        util::ObjectDataLabel( "1DA(omega,alpha)"),
        false,
        "Sequence order should be considered"
      );
      BCL_ExampleIndirectCheck
      (
        util::ObjectDataLabel( "1DA(alpha=5,omega=10)"),
        util::ObjectDataLabel( "1DA(omega=10,alpha=5)"),
        "Key order in maps is irrelevant"
      );
      BCL_ExampleCheck( quoted_stuff.GetNumberArguments(), 0);
      BCL_ExampleCheck( quoted_stuff.GetValue(), "stuff(2)");

      // test the ability to parse a more complex label
      util::ObjectDataLabel data_label( "ResilientPropagation( hidden architecture(16),objective function=RMSD)");
      if( BCL_ExampleCheck( data_label.GetNumberArguments(), 2))
      {
        // check that the primary object itself was not labelled
        BCL_ExampleCheck( data_label.GetName(), "");
        // check the primary object string
        BCL_ExampleCheck( data_label.GetValue(), "ResilientPropagation");

        // check the arguments
        BCL_ExampleCheck( data_label.GetArgument( 0).GetName(), "");
        BCL_ExampleCheck( data_label.GetArgument( 0).GetValue(), "hidden architecture");
        BCL_ExampleCheck( data_label.GetArgument( 0).GetNumberArguments(), 1);
        BCL_ExampleCheck( data_label.GetArgument( 0).GetArgument( 0).ToString(), "16");
        BCL_ExampleCheck( data_label.GetArgument( 1).GetName(), "objective function");
        BCL_ExampleCheck( data_label.GetArgument( 1).GetValue(), "RMSD");
        BCL_ExampleCheck( data_label.GetArgument( 1).GetNumberArguments(), 0);
      }

      // check names with delimiters
      BCL_ExampleCheck
      (
        util::ObjectDataLabel( "\"C=C-O-H,B\"=\"BondType=Aromatic\"(Single)").ToNamedString(),
        "\"C=C-O-H,B\"=\"BondType=Aromatic\"(Single)"
      );

      BCL_ExampleCheck( name_with_delimiter.GetName(), "C=O");
      BCL_ExampleCheck( name_with_delimiter.GetValue(), "BondType=Aromatic");

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

  }; //end ExampleUtilObjectDataLabel

  const ExampleClass::EnumType ExampleUtilObjectDataLabel::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilObjectDataLabel())
  );

} // namespace bcl
