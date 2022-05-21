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
#include "assemble/bcl_assemble_collector_all_possible_domains.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_all_possible_domains.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorAllPossibleDomains :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorAllPossibleDomains *Clone() const
    {
      return new ExampleAssembleCollectorAllPossibleDomains( *this);
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
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test the default constructor
      BCL_MessageStd( "testing default constructor");
      assemble::CollectorAllPossibleDomains def_construct;

      // test the constructor from a domain size and bool to preserve protein model domains
      BCL_MessageStd( "testing constructor to force model connectivity");
      assemble::CollectorAllPossibleDomains complete_construct( math::Range< size_t>( 2, 2), true);

      // test the clone function
      BCL_MessageStd( "testing Clone function");
      util::ShPtr< assemble::CollectorAllPossibleDomains> clone_construct( def_construct.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< assemble::CollectorAllPossibleDomains>(), clone_construct->GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

      // test Collect function using default domain size
      const size_t expected_nr_domains( 5);
      util::ShPtrVector< assemble::Domain> domains( complete_construct.Collect( protein_model));
      BCL_ExampleIndirectCheck( domains.GetSize(), expected_nr_domains, "Collect");

      // locate the SSEs that should be in the first domain
      const assemble::SSE& sse_40_45( *assemble::LocatorSSE( 'A', 40, 45).Locate( protein_model));
      const assemble::SSE& sse_64_72( *assemble::LocatorSSE( 'A', 64, 72).Locate( protein_model));

      // check that they are in the first domain
      BCL_Example_Check
      (
        domains.FirstElement()->GetNumberSSEs() == 2 &&
        domains.FirstElement()->DoesContain( sse_40_45) &&
        domains.FirstElement()->DoesContain( sse_64_72),
        "The first domain has wrong SSEs, it should have 40-45 and 64-72, but instead it is\n" +
          util::Format()( domains.FirstElement()->GetIdentification())
      );

      // test Collect function preserving protein model domains
      const size_t expected_nr_domains_b( 5);
      util::ShPtrVector< assemble::Domain> domains_b( complete_construct.Collect( protein_model));
      BCL_ExampleIndirectCheck( domains_b.GetSize(), expected_nr_domains_b, "Collect");

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( def_construct);

      // read the object back in
      assemble::CollectorAllPossibleDomains collector_read;
      ReadBCLObject( collector_read);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorAllPossibleDomains

  const ExampleClass::EnumType ExampleAssembleCollectorAllPossibleDomains::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorAllPossibleDomains())
  );
  
} // namespace bcl
