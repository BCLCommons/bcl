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
#include "assemble/bcl_assemble_sse_pool_agreement.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_pool_agreement.cpp
  //!
  //! @author woetzen, karakam
  //! @date Jul 18, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEPoolAgreement :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEPoolAgreement *Clone() const
    {
      return new ExampleAssembleSSEPoolAgreement( *this);
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

      // read in pdb for 1ubi
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      //build models from pdb sequences
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel native_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // construct native pool from protein model
      const assemble::SSEPool native_pool( native_model.GetSSEs());

      // initialize read
      io::IFStream read;

      // now read the pools for evaluation
      BCL_ExampleMustOpenInputFile
      (
        read, AddExampleInputPathToFilename( e_Biology, "pool_factory_highest_JUFO.bcl")
      );
      assemble::SSEPool pool_a;
      pool_a.ReadSSEPool( read, native_model, 0, 0);
      io::File::CloseClearFStream( read);

      // read in a second pool
      BCL_ExampleMustOpenInputFile
      (
        read, AddExampleInputPathToFilename( e_Biology, "pool_factory_threshold_JUFO.bcl")
      );
      assemble::SSEPool pool_b;
      pool_b.ReadSSEPool( read, native_model, 0, 0);
      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      assemble::SSEPoolAgreement agreement;

      // constructor from a scheme
      assemble::SSEPoolAgreement agreement_b( true, "test");

    /////////////////
    // data access //
    /////////////////

      // check the default scheme
      BCL_ExampleCheck( agreement.GetDefaultScheme(), "sse_pool_agreement");

      // check the scheme
      BCL_ExampleCheck( agreement_b.GetScheme(), "test");

    ////////////////
    // operations //
    ////////////////

      // expected q3
      const double q3_expected_native_vs_a( 74.5098);
      const double q3_expected_native_vs_b( 77.551);
      const double q3_expected_a_vs_b( 95.4545);

      // calculate actual q3
      const double q3_native_vs_a( agreement.Q3Score( native_pool, pool_a));
      BCL_MessageStd( "Native vs A: " + util::Format()( q3_native_vs_a));
      const double q3_a_vs_native( agreement.Q3Score( pool_a, native_pool));
      BCL_MessageStd( "A vs native: " + util::Format()( q3_a_vs_native));
      const double q3_native_vs_b( agreement.Q3Score( native_pool, pool_b));
      BCL_MessageStd( "Native vs B: " + util::Format()( q3_native_vs_b));
      const double q3_b_vs_native( agreement.Q3Score( pool_b, native_pool));
      BCL_MessageStd( "B vs native: " + util::Format()( q3_b_vs_native));
      const double q3_a_vs_b( agreement.Q3Score( pool_a, pool_b));
      BCL_MessageStd( "A vs B: " + util::Format()( q3_a_vs_b));
      const double q3_b_vs_a( agreement.Q3Score( pool_b, pool_a));
      BCL_MessageStd( "B vs A: " + util::Format()( q3_b_vs_a));

      // compare the agreement scores calculated and the expected scores
      BCL_ExampleCheckWithinTolerance( q3_expected_native_vs_a, q3_native_vs_a, 0.001);
      BCL_ExampleCheckWithinTolerance( q3_expected_native_vs_a, q3_a_vs_native, 0.001);
      BCL_ExampleCheckWithinTolerance( q3_expected_native_vs_b, q3_native_vs_b, 0.001);
      BCL_ExampleCheckWithinTolerance( q3_expected_native_vs_b, q3_b_vs_native, 0.001);
      BCL_ExampleCheckWithinTolerance( q3_expected_a_vs_b, q3_a_vs_b, 0.001);
      BCL_ExampleCheckWithinTolerance( q3_expected_a_vs_b, q3_b_vs_a, 0.001);

      // calculate and test the values of different agreement scores
      const double expected_native_vs_a( 10.0858);
      const double expected_a_vs_native( 5.2575);
      const double expected_native_vs_b( 11.9013);
      const double expected_b_vs_native( 5.66296);
      const double expected_a_vs_b( 6.76157);
      const double expected_b_vs_a( 3.72736);

      // calculate the actual agreement scores
      const double agreement_native_vs_a( agreement.AgreementToTemplate( native_pool, pool_a));
      BCL_MessageStd( "Native vs A: " + util::Format()( agreement_native_vs_a));
      const double agreement_a_vs_native( agreement.AgreementToTemplate( pool_a, native_pool));
      BCL_MessageStd( "A vs native: " + util::Format()( agreement_a_vs_native));
      const double agreement_native_vs_b( agreement.AgreementToTemplate( native_pool, pool_b));
      BCL_MessageStd( "Native vs B: " + util::Format()( agreement_native_vs_b));
      const double agreement_b_vs_native( agreement.AgreementToTemplate( pool_b, native_pool));
      BCL_MessageStd( "B vs native: " + util::Format()( agreement_b_vs_native));
      const double agreement_a_vs_b( agreement.AgreementToTemplate( pool_a, pool_b));
      BCL_MessageStd( "A vs B: " + util::Format()( agreement_a_vs_b));
      const double agreement_b_vs_a( agreement.AgreementToTemplate( pool_b, pool_a));
      BCL_MessageStd( "B vs A: " + util::Format()( agreement_b_vs_a));

      // compare the agreement scores calculated and the expected scores
      BCL_ExampleCheckWithinTolerance( expected_native_vs_a, agreement_native_vs_a, 0.001);
      BCL_ExampleCheckWithinTolerance( expected_a_vs_native, agreement_a_vs_native, 0.001);
      BCL_ExampleCheckWithinTolerance( expected_native_vs_b, agreement_native_vs_b, 0.001);
      BCL_ExampleCheckWithinTolerance( expected_b_vs_native, agreement_b_vs_native, 0.001);
      BCL_ExampleCheckWithinTolerance( expected_a_vs_b, agreement_a_vs_b, 0.001);
      BCL_ExampleCheckWithinTolerance( expected_b_vs_a, agreement_b_vs_a, 0.001);

    ///////////////
    // operators //
    ///////////////

      // expected agreement by operator, which is symmetric
      const double expected_ab( expected_a_vs_b + expected_b_vs_a);

      // calculate and check
      const double agreement_ab( agreement( pool_a, pool_b));
      BCL_ExampleCheckWithinTolerance( expected_ab, agreement_ab, 0.001);
      BCL_MessageStd( "AB agreement: " + util::Format()( agreement_ab));

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( agreement_b);
      assemble::SSEPoolAgreement agreement_read;
      ReadBCLObject( agreement_read);
      BCL_ExampleCheck( agreement_read.GetScheme(), agreement_b.GetScheme());
      BCL_ExampleCheck( agreement_read.AgreementToTemplate( pool_a, pool_b), agreement_b.AgreementToTemplate( pool_a, pool_b));

    //////////////////////
    // helper functions //
    //////////////////////

      // check the function overlap for single SSEs
      const util::SiPtr< const assemble::SSE> helix_23_34( assemble::LocatorSSE( 'A', 23, 34).Locate( native_model));
      const util::SiPtr< const assemble::SSE> strand_10_17( assemble::LocatorSSE( 'A', 10, 17).Locate( native_model));
      const util::SiPtr< const assemble::SSE> strand_12_17( assemble::LocatorSSE( 'A', 12, 17).Locate( pool_a));

      // initialize expected overlaps
      const storage::VectorND< 2, int> expected_overlap_a( 13, -17);
      const storage::VectorND< 2, int> expected_overlap_b( -13, 17);
      const storage::VectorND< 2, int> expected_overlap_c( 2, 0);
      const storage::VectorND< 2, int> expected_overlap_d( -2, 0);

      // calculate overlaps
      const storage::VectorND< 2, int> overlap_a( assemble::SSEPoolAgreement::Overlap( *strand_10_17, *helix_23_34));
      BCL_MessageVrb( "overlap_a: " + util::Format()( overlap_a));
      const storage::VectorND< 2, int> overlap_b( assemble::SSEPoolAgreement::Overlap( *helix_23_34, *strand_10_17));
      BCL_MessageVrb( "overlap_b: " + util::Format()( overlap_b));
      const storage::VectorND< 2, int> overlap_c( assemble::SSEPoolAgreement::Overlap( *strand_10_17, *strand_12_17));
      BCL_MessageVrb( "overlap_c: " + util::Format()( overlap_c));
      const storage::VectorND< 2, int> overlap_d( assemble::SSEPoolAgreement::Overlap( *strand_12_17, *strand_10_17));
      BCL_MessageVrb( "overlap_d: " + util::Format()( overlap_d));

      // check the overlaps
      BCL_ExampleCheck( overlap_a, expected_overlap_a);
      BCL_ExampleCheck( overlap_b, expected_overlap_b);
      BCL_ExampleCheck( overlap_c, expected_overlap_c);
      BCL_ExampleCheck( overlap_d, expected_overlap_d);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEPoolAgreement

  const ExampleClass::EnumType ExampleAssembleSSEPoolAgreement::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEPoolAgreement())
  );

} // namespace bcl
