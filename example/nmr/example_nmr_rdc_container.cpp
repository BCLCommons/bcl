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
#include "nmr/bcl_nmr_rdc_container.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_nmr_rdc_container.cpp
  //!
  //! @author weinerbe
  //! @date Mar 23, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleNmrRDCContainer :
    public ExampleInterface
  {
  public:

    ExampleNmrRDCContainer *Clone() const
    {
      return new ExampleNmrRDCContainer( *this);
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
      nmr::RDCContainer def_construct;
      BCL_ExampleIndirectCheck( def_construct.GetExperimentalValues().IsEmpty(), true, "default constructor");

      // test constructor from values
      const storage::Vector< double> exp_values( storage::Vector< double>::Create( 1.5, -0.4, 2.8));
      const storage::Vector< double> calc_values( storage::Vector< double>::Create( 1.2, -0.5, 2.8));
      const nmr::RDCContainer value_construct( exp_values, calc_values);

    /////////////////
    // data access //
    /////////////////

      // test GetExperimentalValues
      BCL_ExampleCheck( value_construct.GetExperimentalValues(), exp_values);

      // test GetCalculatedlValues
      BCL_ExampleCheck( value_construct.GetCalculatedlValues(), calc_values);

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( value_construct);
      ReadBCLObject( def_construct);
      BCL_ExampleIndirectCheck
      (
        def_construct.GetExperimentalValues() == exp_values && def_construct.GetCalculatedlValues() == calc_values,
        true,
        "read and write"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleNmrRDCContainer

  const ExampleClass::EnumType ExampleNmrRDCContainer::s_Instance
  (
    GetExamples().AddEnum( ExampleNmrRDCContainer())
  );
  
} // namespace bcl
