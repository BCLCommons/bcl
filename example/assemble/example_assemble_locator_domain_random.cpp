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
#include "assemble/bcl_assemble_locator_domain_random.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_locator_domain_random.cpp
  //!
  //! @author karakam
  //! @date Feb 3, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleLocatorDomainRandom :
    public ExampleInterface
  {
  public:

    ExampleAssembleLocatorDomainRandom *Clone() const
    {
      return new ExampleAssembleLocatorDomainRandom( *this);
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
      // pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1J27A.pdb"));

      // initialize min sse sizes
      storage::Map< biol::SSType, size_t> sse_sizes;
      sse_sizes[ biol::GetSSTypes().HELIX] = 5;
      sse_sizes[ biol::GetSSTypes().STRAND] = 5;
      sse_sizes[ biol::GetSSTypes().COIL] = 999;

      // read protein model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, sse_sizes));

      // ranges
      math::Range< size_t> range_default( 1, 1);
      math::Range< size_t> range_a( 2, 2);
      math::Range< size_t> range_b( 2, 4);
      math::Range< size_t> range_c( 20, 30);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      assemble::LocatorDomainRandom locator_default;
      BCL_ExampleCheck( locator_default.GetDomainSizeRange(), range_default);
      BCL_ExampleCheck( locator_default.GetSSType(), biol::GetSSTypes().COIL);

      // constructor from a range and sstype
      assemble::LocatorDomainRandom locator_a( range_a, biol::GetSSTypes().HELIX);
      assemble::LocatorDomainRandom locator_b( range_b, biol::GetSSTypes().STRAND);
      BCL_ExampleCheck( locator_b.GetDomainSizeRange(), range_b);
      BCL_ExampleCheck( locator_b.GetSSType(), biol::GetSSTypes().STRAND);
      assemble::LocatorDomainRandom locator_c( range_c, biol::GetSSTypes().STRAND);
      BCL_ExampleCheck( locator_c.GetDomainSizeRange(), range_c);
      BCL_ExampleCheck( locator_c.GetSSType(), biol::GetSSTypes().STRAND);

    /////////////////
    // data access //
    /////////////////

      // test GetDomainSizeRange()
      BCL_MessageStd( "Test GetDomainSizeRange()");
      BCL_ExampleCheck( locator_a.GetDomainSizeRange(), range_a);

      // test GetSSType()
      BCL_MessageStd( "Test GetSSType()");
      BCL_ExampleCheck( locator_a.GetSSType(), biol::GetSSTypes().HELIX);

    ////////////////
    // operations //
    ////////////////

      // test Locate function
      BCL_MessageStd( "Test Locate function locator_a");
      util::ShPtr< assemble::Domain> sp_domain_a( locator_a.Locate( model));
      BCL_ExampleCheck( sp_domain_a.IsDefined(), true);
      BCL_ExampleCheck( sp_domain_a->GetNumberSSEs(), 2);
      BCL_ExampleCheck( sp_domain_a->GetNumberSSE( biol::GetSSTypes().HELIX), 2);

      // test locate function
      BCL_MessageStd( "Test Locate function with locator_b");
      util::ShPtr< assemble::Domain> sp_domain_b( locator_b.Locate( model));
      BCL_ExampleCheck( sp_domain_b.IsDefined(), true);
      BCL_ExampleCheck
      (
        sp_domain_b->GetNumberSSEs(), sp_domain_b->GetNumberSSE( biol::GetSSTypes().STRAND)
      );
      BCL_ExampleCheck
      (
        range_b.IsWithin( sp_domain_b->GetNumberSSE( biol::GetSSTypes().STRAND)), true
      );

      // test locate function
      BCL_MessageStd( "Test Locate function with locator_c");
      util::ShPtr< assemble::Domain> sp_domain_c( locator_c.Locate( model));
      BCL_ExampleCheck( sp_domain_c.IsDefined(), false);

    //////////////////////
    // input and output //
    //////////////////////

      // test input output
      WriteBCLObject( locator_a);
      // initialize new object and read
      assemble::LocatorDomainRandom locator_read;
      ReadBCLObject( locator_read);
      // compare
      BCL_ExampleIndirectCheck
      (
        locator_read.GetDomainSizeRange() == locator_a.GetDomainSizeRange() &&
        locator_read.GetSSType() == locator_a.GetSSType(),
        true,
        "The read function failed"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleLocatorDomainRandom

  const ExampleClass::EnumType ExampleAssembleLocatorDomainRandom::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleLocatorDomainRandom())
  );

} // namespace bcl
