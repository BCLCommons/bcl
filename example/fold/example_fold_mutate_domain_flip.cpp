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
#include "fold/bcl_fold_mutate_domain_flip.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_domain_flip.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateDomainFlip :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateDomainFlip *Clone() const
    {
      return new ExampleFoldMutateDomainFlip( *this);
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

      // construct min_sse_sizes
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 6;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 999;

      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2yv8_ideal.pdb"));

      // get the protein model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes));

      // join the adjacent strand if any exists
      model.Join( biol::GetSSTypes().STRAND, false);

      // create a pointer on chain
      util::ShPtr< assemble::Chain> chain( model.GetChain( 'A'));

      // collect the sheets and create pointers for conveniance
      util::ShPtrVector< assemble::Domain> sheets( assemble::CollectorSheet().Collect( model));
      util::ShPtr< assemble::Domain> sheet_a( sheets( 1));
      util::ShPtr< assemble::Domain> sheet_b( sheets( 0));

      // print out the contents of the sheets
      BCL_MessageStd( "Sheet#1:\n" + sheet_a->GetTopology()->GetOrderedIdentification());
      BCL_MessageStd( "Sheet#2:\n" + sheet_b->GetTopology()->GetOrderedIdentification());

      // construct set of axes to flip
      storage::Set< coord::Axis> xyz_axes;
      xyz_axes.Insert( coord::GetAxes().e_X);
      xyz_axes.Insert( coord::GetAxes().e_Y);
      xyz_axes.Insert( coord::GetAxes().e_Z);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      BCL_MessageStd( "test default constructor");
      fold::MutateDomainFlip mutate_def;

      BCL_MessageStd( "test constructor from an axis, external");
      fold::MutateDomainFlip mutate_z_external( coord::GetAxes().e_Z, false);

      BCL_MessageStd( "test constructor from a set of axes, internal, subset");
      fold::MutateDomainFlip mutate_xyz_internal_subset( xyz_axes, true, false, true);

      BCL_MessageStd( "test constructor from an axis, internal, all");
      fold::MutateDomainFlip mutate_x_internal( coord::GetAxes().e_X, true, true, false);

      BCL_MessageStd( "test clone constructor");
      util::ShPtr< fold::MutateDomainFlip> sp_mutate( mutate_z_external.Clone());
      BCL_Example_Check
      (
        sp_mutate->GetFlipAxes().InternalData() == mutate_z_external.GetFlipAxes().InternalData() &&
        sp_mutate->GetFlipInternal() == mutate_z_external.GetFlipInternal() &&
        sp_mutate->GetFlipAll() == mutate_z_external.GetFlipAll() &&
        sp_mutate->GetUseDifferentFlipAxes() == mutate_z_external.GetUseDifferentFlipAxes(),
        "The cloned mutate is different!\n" + util::Format()( *sp_mutate) + "\nvs\n" + util::Format()( mutate_z_external)
      );

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_MessageStd
      (
        "This class has the following identifier " + mutate_def.GetClassIdentifier()
      );

      // test GetFlipAxes
      BCL_MessageStd( "test GetFlipAxes()");
      BCL_Example_Check
      (
        mutate_xyz_internal_subset.GetFlipAxes().InternalData() == xyz_axes.InternalData(),
        "GetFlipAxes() should return " + util::Format()( xyz_axes) +
        " not " + util::Format()( mutate_xyz_internal_subset.GetFlipAxes())
      );

      // test GetFlipInternal
      BCL_MessageStd( "test GetFlipInternal()");
      BCL_Example_Check
      (
        mutate_xyz_internal_subset.GetFlipInternal(),
        "GetFlipInternal() should return true not false!"
      );

      // test GetFlipAll
      BCL_MessageStd( "test GetFlipAll()");
      BCL_Example_Check
      (
        !mutate_xyz_internal_subset.GetFlipAll(),
        "GetFlipAll() should return false not true!"
      );

      // test GetUseDifferentFlipAxes
      BCL_MessageStd( "test GetUseDifferentFlipAxes()");
      BCL_Example_Check
      (
        mutate_xyz_internal_subset.GetUseDifferentFlipAxes(),
        "GetUseDifferentFlipAxes() should be return true not false!!"
      );

    ///////////////
    // operators //
    ///////////////

    ///////////////////////////////
    // mutate_z_external sheet_a //
    ///////////////////////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_z_external first sheet");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_a( mutate_z_external( *sheet_a));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_a.GetArgument().IsDefined(), true, "The mutate returned an empty argument!!"
      );

      // construct a chain
      assemble::Chain chain_a( chain->GetSequence(), *mutate_result_a.GetArgument());

      // write out the chain
      Proteins::WriteChainToPDB( chain_a, AddExampleOutputPathToFilename( mutate_def, "mutate_domain_flip_a.pdb"));

    ///////////////////////////////
    // mutate_z_external sheet_b //
    ///////////////////////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_z_external second sheet");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_b( mutate_z_external( *sheet_b));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_b.GetArgument().IsDefined(), true, "The mutate returned an empty argument!!"
      );

      // construct a chain
      assemble::Chain chain_b( chain->GetSequence(), *mutate_result_b.GetArgument());

      // write out the chain
      Proteins::WriteChainToPDB( chain_b, AddExampleOutputPathToFilename( mutate_def, "mutate_domain_flip_b.pdb"));

    ///////////////////////////////
    // mutate_x_internal sheet_a //
    ///////////////////////////////

      // test the operator()
      BCL_MessageStd( "testing operator() with mutate_x_internal");

      // test with a single mutate
      math::MutateResult< assemble::Domain> mutate_result_c( mutate_x_internal( *sheet_a));

      // make sure the return result is not empty
      BCL_ExampleIndirectAssert
      (
        mutate_result_c.GetArgument().IsDefined(), true, "The mutate returned an empty argument!!"
      );

       // construct a chain
       assemble::Chain chain_c( chain->GetSequence(), *mutate_result_c.GetArgument());

       // write out the model
       Proteins::WriteChainToPDB( chain_c, AddExampleOutputPathToFilename( mutate_def, "mutate_domain_flip_c.pdb"));

    ////////////////////////////////////////
    // mutate_xyz_internal_subset sheet_a //
    ////////////////////////////////////////

       // test the operator()
       BCL_MessageStd( "testing operator() with mutate_xyz_internal_subset");

       // test with a single mutate
       math::MutateResult< assemble::Domain> mutate_result_d( mutate_xyz_internal_subset( *sheet_a));

       // make sure the return result is not empty
       BCL_ExampleIndirectAssert
       (
         mutate_result_d.GetArgument().IsDefined(), true, "The mutate returned an empty argument!!"
       );

       // construct a chain
       assemble::Chain chain_d( chain->GetSequence(), *mutate_result_d.GetArgument());

       // write out the model
       Proteins::WriteChainToPDB( chain_d, AddExampleOutputPathToFilename( mutate_def, "mutate_domain_flip_d.pdb"));

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "testing read write")
      WriteBCLObject( mutate_z_external);
      fold::MutateDomainFlip mutate_read;
      ReadBCLObject( mutate_read);
      BCL_Example_Check
      (
        mutate_read.GetFlipAxes().InternalData() == mutate_z_external.GetFlipAxes().InternalData() &&
        mutate_read.GetFlipInternal() == mutate_z_external.GetFlipInternal() &&
        mutate_read.GetFlipAll() == mutate_z_external.GetFlipAll() &&
        mutate_read.GetUseDifferentFlipAxes() == mutate_z_external.GetUseDifferentFlipAxes(),
        "The read mutate is different\n" + util::Format()( mutate_read) + "\nvs\n" + util::Format()( mutate_z_external)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateDomainFlip

  const ExampleClass::EnumType ExampleFoldMutateDomainFlip::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateDomainFlip())
  );

} // namespace bcl
