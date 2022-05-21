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
#include "fold/bcl_fold_handler_locator_loop_segment.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_locator_loop_segment.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_handler_locator_loop_segment.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldHandlerLocatorLoopSegment :
    public ExampleInterface
  {
  public:

    ExampleFoldHandlerLocatorLoopSegment *Clone() const
    {
      return new ExampleFoldHandlerLocatorLoopSegment( *this);
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

      {
        const std::string loop_segment_filename( AddExampleInputPathToFilename( e_Fold, "handler_loop_segment_rigid.txt"));
        fold::HandlerLocatorLoopSegment handler;
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, loop_segment_filename);
        fold::LocatorLoopSegment created_locator( handler.HandleRead( read));
        BCL_ExampleCheck
        (
          ( created_locator.GetLocatorSSE().GetSSEID() == storage::VectorND< 2, int>( 35, 38)), true
        );
        BCL_ExampleCheck
        (
          created_locator.GetLocatorSSE().GetChainID() == 'A', true
        );
        BCL_ExampleCheck
        (
          created_locator.IsRigid(), true
        );
      }

      {
        const std::string loop_segment_filename( AddExampleInputPathToFilename( e_Fold, "handler_loop_segment_not_rigid.txt"));
        fold::HandlerLocatorLoopSegment handler;
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, loop_segment_filename);
        fold::LocatorLoopSegment created_locator( handler.HandleRead( read));
        BCL_ExampleCheck
        (
          ( created_locator.GetLocatorSSE().GetSSEID() == storage::VectorND< 2, int>( 35, 38)), true
        );
        BCL_ExampleCheck
        (
          created_locator.GetLocatorSSE().GetChainID() == 'A', true
        );
        BCL_ExampleCheck
        (
          created_locator.IsRigid(), false
        );
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

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

  }; //end ExampleFoldHandlerLocatorLoopSegment

  const ExampleClass::EnumType ExampleFoldHandlerLocatorLoopSegment::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldHandlerLocatorLoopSegment())
  );

} // namespace bcl
