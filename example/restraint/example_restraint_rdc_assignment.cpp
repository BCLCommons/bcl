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
#include "restraint/bcl_restraint_rdc_assignment.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_rdc_assignment.cpp
  //!
  //! @author weinerbe
  //! @date Mar 21, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintRDCAssignment :
    public ExampleInterface
  {
  public:

    ExampleRestraintRDCAssignment *Clone() const
    {
      return new ExampleRestraintRDCAssignment( *this);
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

      // test default constructor
      restraint::RDCAssignment def_construct;

    /////////////////
    // data access //
    /////////////////

      // test GetData
      BCL_ExampleIndirectCheck( def_construct.GetData().GetSize(), 0, "GetData");

    ////////////////
    // operations //
    ////////////////

      // test PushBack
      def_construct.PushBack( linal::Vector3D(), linal::Vector3D(), 1.0);
      def_construct.PushBack( linal::Vector3D(), linal::Vector3D(), 2.0);
      def_construct.PushBack( linal::Vector3D(), linal::Vector3D(), 3.0);
      BCL_ExampleIndirectCheck( def_construct.GetData().FirstElement().Third(), 1.0, "PushBack");

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( def_construct);
      restraint::RDCAssignment read_construct;
      ReadBCLObject( read_construct);
      BCL_ExampleIndirectCheck
      (
        def_construct.GetData().GetSize(),
        read_construct.GetData().GetSize(),
        "read and write"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintRDCAssignment

  const ExampleClass::EnumType ExampleRestraintRDCAssignment::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintRDCAssignment())
  );
  
} // namespace bcl
