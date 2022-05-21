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
#include "fold/bcl_fold_mutate_sse_bend_template.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_compare.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_sse_bend_template.cpp
  //!
  //! @author weinerbe
  //! @date Jul 6, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSSEBendTemplate :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSSEBendTemplate *Clone() const
    {
      return new ExampleFoldMutateSSEBendTemplate( *this);
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
      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));

      // get the SSEs
      util::ShPtr< assemble::SSE>
        sp_helix( Proteins::GetSSE( pdb_filename, 'A', 23, 34, biol::GetAAClasses().e_AABackBone)),
        sp_strand( Proteins::GetSSE( pdb_filename, 'A', 1, 7, biol::GetAAClasses().e_AABackBone));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      fold::MutateSSEBendTemplate def_construct;
      BCL_ExampleCheck( def_construct.GetScheme(), GetStaticClassName< fold::MutateSSEBendTemplate>());

      // test construct from comparison
      const fold::MutateSSEBendTemplate comparison_construct
      (
        assemble::SSEGeometryWithinSizeTolerance( 0, 0),
        "size_0"
      );

    /////////////////
    // data access //
    /////////////////

      // test GetScheme
      BCL_ExampleCheck( comparison_construct.GetScheme(), "size_0");

    ///////////////
    // operators //
    ///////////////

      // test () operator w/ helix
      BCL_MessageStd( "test operator() with helix");
      math::MutateResult< assemble::SSE> result_helix( comparison_construct( *sp_helix));
      BCL_ExampleCheck( result_helix.GetArgument().IsDefined(), true);
      BCL_ExampleCheck( result_helix.GetArgument()->HasDefinedCoordinates(), true);

      BCL_MessageStd( "test operator() with strand");
      math::MutateResult< assemble::SSE> result_strand( comparison_construct( *sp_strand));
      BCL_ExampleCheck( result_strand.GetArgument().IsDefined(), true);
      BCL_ExampleCheck( result_strand.GetArgument()->HasDefinedCoordinates(), true);

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( comparison_construct);
      ReadBCLObject( def_construct);
      BCL_ExampleIndirectCheck( def_construct.GetScheme(), comparison_construct.GetScheme(), "read and write");

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSSEBendTemplate

  const ExampleClass::EnumType ExampleFoldMutateSSEBendTemplate::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSSEBendTemplate())
  );
  
} // namespace bcl
