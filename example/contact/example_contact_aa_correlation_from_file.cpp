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
#include "contact/bcl_contact_aa_correlation_from_file.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_aa_correlation_from_file.cpp
  //!
  //! @author teixeipl, mendenjl
  //! @date Feb 11, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactAACorrelationFromFile :
    public ExampleInterface
  {
  public:

    ExampleContactAACorrelationFromFile *Clone() const
    {
      return new ExampleContactAACorrelationFromFile( *this);
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
      // form an AA
      biol::AA amino_acid( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 1, 1)));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      contact::AACorrelationFromFile aa_info;

      // clone
      util::ShPtr< contact::AACorrelationFromFile> sp_aa_info( aa_info.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetLength function
      BCL_ExampleCheck( aa_info.GetSizeOfFeatures(), 1);

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

  }; //end ExampleContactAACorrelationFromFile

  const ExampleClass::EnumType ExampleContactAACorrelationFromFile::s_Instance
  (
    GetExamples().AddEnum( ExampleContactAACorrelationFromFile())
  );

} // namespace bcl
