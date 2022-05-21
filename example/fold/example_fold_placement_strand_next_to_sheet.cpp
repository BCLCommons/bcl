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
#include "fold/bcl_fold_placement_strand_next_to_sheet.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_placement_strand_next_to_sheet.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldPlacementStrandNextToSheet :
    public ExampleInterface
  {
  public:

    ExampleFoldPlacementStrandNextToSheet *Clone() const
    {
      return new ExampleFoldPlacementStrandNextToSheet( *this);
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

      // initialize pdb filename for 1ubi
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2yv8_ideal.pdb"));
      // initialize min sses sizes
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 5;

      // get the model
      assemble::ProteinModel original_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );

      // create a handle on the chain
      util::ShPtr< assemble::Chain> sp_orig_chain( original_model.GetChains().FirstElement());
      // create a handle on the sequence
      util::ShPtr< biol::AASequence> sp_orig_sequence( sp_orig_chain->GetSequence());

      // create locators for the strands
      assemble::LocatorSSE locator_38_45( 'A', 38, 45);
      assemble::LocatorSSE locator_53_59( 'A', 53, 59);
      assemble::LocatorSSE locator_68_76( 'A', 68, 76);
      assemble::LocatorSSE locator_82_89( 'A', 82, 89);
      assemble::LocatorSSE locator_109_116( 'A', 109, 116);
      assemble::LocatorSSE locator_120_125( 'A', 120, 125);
      assemble::LocatorSSE locator_128_134( 'A', 128, 134);
      assemble::LocatorSSE locator_144_149( 'A', 144, 149);
      assemble::LocatorSSE locator_152_159( 'A', 152, 159);

      // get ShPtrs to the original SSEs
      util::SiPtr< const assemble::SSE> strand_38_45( locator_38_45.Locate( original_model));
      util::SiPtr< const assemble::SSE> strand_53_59( locator_53_59.Locate( original_model));
      util::SiPtr< const assemble::SSE> strand_68_76( locator_68_76.Locate( original_model));
      util::SiPtr< const assemble::SSE> strand_82_89( locator_82_89.Locate( original_model));
      util::SiPtr< const assemble::SSE> strand_109_116( locator_109_116.Locate( original_model));
      util::SiPtr< const assemble::SSE> strand_120_125( locator_120_125.Locate( original_model));
      util::SiPtr< const assemble::SSE> strand_128_134( locator_120_125.Locate( original_model));
      util::SiPtr< const assemble::SSE> strand_144_149( locator_144_149.Locate( original_model));
      util::SiPtr< const assemble::SSE> strand_152_159( locator_152_159.Locate( original_model));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      fold::PlacementStrandNextToSheet placement;

      // test constructor with a flip probability
      fold::PlacementStrandNextToSheet placement_no_flip( 0.0);
      fold::PlacementStrandNextToSheet placement_always_flip( 1.0);

      // test copy constructor
      fold::PlacementStrandNextToSheet placement_copy( placement);

      // test clone constructor
      util::ShPtr< fold::PlacementStrandNextToSheet> placement_clone( placement.Clone());

    /////////////////
    // data access //
    /////////////////

      // example check
      BCL_ExampleCheck( placement.GetClassIdentifier(), GetStaticClassName( placement));

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////////
    // test with one strand //
    //////////////////////////
      {
        BCL_MessageStd( "========>  test with one strand");

      /////////////
      // no flip //
      /////////////

        // test with one_strand
        // create a model that has only one strand
        assemble::ProteinModel model_one_strand;
        model_one_strand.Insert( util::ShPtr< assemble::Chain>( new assemble::Chain( sp_orig_sequence)));
        model_one_strand.Insert( util::ShPtr< assemble::SSE>( strand_109_116->HardCopy()));
        // make a copy of the strand to be placed
        util::ShPtr< assemble::SSE> strand_38_45_copy_a( strand_38_45->HardCopy());
        strand_38_45_copy_a->Transform( math::Inverse( strand_38_45_copy_a->GetOrientation()));

        // write out the pdb for this model
        Proteins::WriteModelToPDB
        (
          model_one_strand,
          AddExampleOutputPathToFilename( placement, "placement_strand_next_to_sheet_one.pdb")
        );

        // call the Place and store the result
        storage::Pair< math::TransformationMatrix3D, bool> transformation_a
        (
          placement_no_flip.Place( *strand_38_45_copy_a, model_one_strand)
        );

        // make sure a placement was found
        BCL_Assert( transformation_a.Second(), "The placement failed for model_one_strand");

        // apply the transformation and add it to model
        strand_38_45_copy_a->Transform( transformation_a.First());
        model_one_strand.Insert( strand_38_45_copy_a);

        // check that the strand was added correctly
        BCL_Example_Check
        (
          model_one_strand.GetNumberSSEs() == 2 && model_one_strand.DoesContain( *strand_38_45_copy_a),
          "The model now should have two strands including strand_38_45 but has " +
          util::Format()( model_one_strand.GetNumberSSEs())
        );

        // write out the pdb for this model
        Proteins::WriteModelToPDB
        (
          model_one_strand,
          AddExampleOutputPathToFilename( placement, "placement_strand_next_to_sheet_one_no_flip.pdb")
        );

      ///////////////
      // with flip //
      ///////////////

        // create a model that has only one strand
        assemble::ProteinModel model_one_strand_copy;
        model_one_strand_copy.Insert( util::ShPtr< assemble::Chain>( new assemble::Chain( sp_orig_sequence)));
        model_one_strand_copy.Insert( util::ShPtr< assemble::SSE>( strand_109_116->HardCopy()));

        // make a copy of the strand to be placed
        util::ShPtr< assemble::SSE> strand_38_45_copy_b( strand_38_45->HardCopy());
        strand_38_45_copy_b->SetToIdealConformationAtOrigin();

        // call the Place and store the result
        transformation_a = placement_always_flip.Place( *strand_38_45_copy_b, model_one_strand_copy);

        // make sure a placement was found
        BCL_Assert( transformation_a.Second(), "The placement failed for model_one_strand");

        // apply the transformation and add it to model
        strand_38_45_copy_b->Transform( transformation_a.First());
        model_one_strand_copy.Insert( strand_38_45_copy_b);

        // check that the strand was added correctly
        BCL_Example_Check
        (
          model_one_strand_copy.GetNumberSSEs() == 2 && model_one_strand_copy.DoesContain( *strand_38_45_copy_b),
          "The model now should have two strands including strand_38_45 but has " +
          util::Format()( model_one_strand_copy.GetNumberSSEs())
        );

        // write out the pdb for this model
        Proteins::WriteModelToPDB
        (
          model_one_strand_copy,
          AddExampleOutputPathToFilename( placement, "placement_strand_next_to_sheet_one_with_flip.pdb")
        );
      }

    ///////////////////////////
    // test with two strands //
    ///////////////////////////
      {
      /////////////
      // no flip //
      /////////////

        BCL_MessageStd( "========>  test with two strands");
        // create a model that has two strands
        assemble::ProteinModel model_two_strands;
        model_two_strands.Insert( util::ShPtr< assemble::Chain>( new assemble::Chain( sp_orig_sequence)));
        model_two_strands.Insert( util::ShPtr< assemble::SSE>( strand_38_45->HardCopy()));
        model_two_strands.Insert( util::ShPtr< assemble::SSE>( strand_109_116->HardCopy()));
        // write out the pdb for this model
        Proteins::WriteModelToPDB
        (
          model_two_strands,
          AddExampleOutputPathToFilename( placement, "placement_strand_next_to_sheet_two.pdb")
        );
        //
        // make a copy of the strand to be placed
        util::ShPtr< assemble::SSE> strand_53_59_copy_a( strand_53_59->HardCopy());
        strand_53_59_copy_a->SetToIdealConformationAtOrigin();

        // call the Place and store the result
        storage::Pair< math::TransformationMatrix3D, bool> transformation_b
        (
          placement_no_flip.Place( *strand_53_59_copy_a, model_two_strands)
        );

        // make sure a placement was found
        BCL_Assert( transformation_b.Second(), "The placement failed for model_two_strands");

        // apply the transformation and add it to model
        strand_53_59_copy_a->Transform( transformation_b.First());
        model_two_strands.Insert( strand_53_59_copy_a);

        // check that the strand was added correctly
        BCL_Example_Check
        (
          model_two_strands.GetNumberSSEs() == 3 && model_two_strands.DoesContain( *strand_53_59_copy_a),
          "The model now should have three strands including strand_53_59 but has " +
          util::Format()( model_two_strands.GetNumberSSEs())
        );

        // write out the pdb for this model
        Proteins::WriteModelToPDB
        (
          model_two_strands,
          AddExampleOutputPathToFilename( placement, "placement_strand_next_to_sheet_two_no_flip.pdb")
        );

      ///////////////
      // with flip //
      ///////////////

        // make a hardcopy of the model
        assemble::ProteinModel model_two_strands_copy;
        model_two_strands_copy.Insert( util::ShPtr< assemble::Chain>( new assemble::Chain( sp_orig_sequence)));
        model_two_strands_copy.Insert( util::ShPtr< assemble::SSE>( strand_38_45->HardCopy()));
        model_two_strands_copy.Insert( util::ShPtr< assemble::SSE>( strand_109_116->HardCopy()));

        // make a copy of the strand to be placed
        util::ShPtr< assemble::SSE> strand_53_59_copy_b( strand_53_59->HardCopy());
        strand_53_59_copy_b->SetToIdealConformationAtOrigin();

        // call the Place and store the result
        transformation_b = placement_always_flip.Place( *strand_53_59_copy_b, model_two_strands_copy);

        // make sure a placement was found
        BCL_Assert( transformation_b.Second(), "The placement failed for model_two_strands");

        // apply the transformation and add it to model
        strand_53_59_copy_b->Transform( transformation_b.First());
        model_two_strands_copy.Insert( strand_53_59_copy_b);

        // check that the strand was added correctly
        BCL_Example_Check
        (
          model_two_strands_copy.GetNumberSSEs() == 3 && model_two_strands_copy.DoesContain( *strand_53_59_copy_b),
          "The model now should have three strands including strand_53_59 but has " +
          util::Format()( model_two_strands_copy.GetNumberSSEs())
        );

        // write out the pdb for this model
        Proteins::WriteModelToPDB
        (
          model_two_strands_copy,
          AddExampleOutputPathToFilename( placement, "placement_strand_next_to_sheet_two_with_flip.pdb")
        );
      }

    /////////////////////////////
    // test with three strands //
    /////////////////////////////
      {
        BCL_MessageStd( "========> Test with three strands");

      /////////////
      // no flip //
      /////////////

        // test with three_strands
        // create a model that has three strands
        assemble::ProteinModel model_three_strands;
        model_three_strands.Insert( util::ShPtr< assemble::Chain>( new assemble::Chain( sp_orig_sequence)));
        model_three_strands.Insert( util::ShPtr< assemble::SSE>( strand_53_59->HardCopy()));
        model_three_strands.Insert( util::ShPtr< assemble::SSE>( strand_68_76->HardCopy()));
        model_three_strands.Insert( util::ShPtr< assemble::SSE>( strand_144_149->HardCopy()));

        // write out the pdb for this model
        Proteins::WriteModelToPDB
        (
          model_three_strands,
          AddExampleOutputPathToFilename( placement, "placement_strand_next_to_sheet_three.pdb")
        );

        // make a copy of the strand to be placed
        util::ShPtr< assemble::SSE> strand_82_89_copy_a( strand_82_89->HardCopy());
        strand_82_89_copy_a->SetToIdealConformationAtOrigin();

        // call the Place and store the result
        storage::Pair< math::TransformationMatrix3D, bool> transformation_c
        (
          placement_no_flip.Place( *strand_82_89_copy_a, model_three_strands)
        );

        // make sure a placement was found
        BCL_Assert( transformation_c.Second(), "The placement failed for model_three_strands");

        // apply the transformation and add it to model
        strand_82_89_copy_a->Transform( transformation_c.First());
        model_three_strands.Insert( strand_82_89_copy_a);

        // check that the strand was added correctly
        BCL_Example_Check
        (
          model_three_strands.GetNumberSSEs() == 4 && model_three_strands.DoesContain( *strand_82_89_copy_a),
          "The model now should have four strands including strand_82_89 but has " +
          util::Format()( model_three_strands.GetNumberSSEs())
        );

        // write out the pdb for this model
        Proteins::WriteModelToPDB
        (
          model_three_strands,
          AddExampleOutputPathToFilename( placement, "placement_strand_next_to_sheet_three_no_flip.pdb")
        );

      ///////////////
      // with flip //
      ///////////////

        // create a model that has three strands
        assemble::ProteinModel model_three_strands_copy;
        model_three_strands_copy.Insert( util::ShPtr< assemble::Chain>( new assemble::Chain( sp_orig_sequence)));
        model_three_strands_copy.Insert( util::ShPtr< assemble::SSE>( strand_53_59->HardCopy()));
        model_three_strands_copy.Insert( util::ShPtr< assemble::SSE>( strand_68_76->HardCopy()));
        model_three_strands_copy.Insert( util::ShPtr< assemble::SSE>( strand_144_149->HardCopy()));

        // make a copy of the strand to be placed
        util::ShPtr< assemble::SSE> strand_82_89_copy_b( strand_82_89->HardCopy());
        strand_82_89_copy_b->SetToIdealConformationAtOrigin();

        // call the Place and store the result
        transformation_c = placement_always_flip.Place( *strand_82_89_copy_b, model_three_strands_copy);

        // make sure a placement was found
        BCL_Assert( transformation_c.Second(), "The placement failed for model_three_strands");

        // apply the transformation and add it to model
        strand_82_89_copy_b->Transform( transformation_c.First());
        model_three_strands_copy.Insert( strand_82_89_copy_b);

        // check that the strand was added correctly
        BCL_Example_Check
        (
          model_three_strands_copy.GetNumberSSEs() == 4 && model_three_strands_copy.DoesContain( *strand_82_89_copy_b),
          "The model now should have four strands including strand_82_89 but has " +
          util::Format()( model_three_strands_copy.GetNumberSSEs())
        );

        // write out the pdb for this model
        Proteins::WriteModelToPDB
        (
          model_three_strands_copy,
          AddExampleOutputPathToFilename( placement, "placement_strand_next_to_sheet_three_with_flip.pdb")
        );
      }
    //////////////////////
    // iterative adding //
    //////////////////////
      {
        // get all the strands
        util::SiPtrVector< const assemble::SSE>
        all_strands( original_model.GetSSEs( biol::GetSSTypes().STRAND).SubSiPtrVector( 0, 5));

        // create a model with the first one
        assemble::ProteinModel iterative_model;
        iterative_model.Insert( util::ShPtr< assemble::Chain>( new assemble::Chain( sp_orig_sequence)));
        // insert the first SSE into model and remove from the strands list
        iterative_model.Insert( util::ShPtr< assemble::SSE>( all_strands.FirstElement()->HardCopy()));
        all_strands.Remove( all_strands.Begin());
        // write the model
        Proteins::WriteModelToPDB
        (
          iterative_model, AddExampleOutputPathToFilename( placement, "placement_strand_next_to_sheet_itr_1.pdb")
        );

        // counter
        size_t strand_ctr( 1);
        // now iterate over the remaining strands
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
          strand_itr( all_strands.Begin()), strand_itr_end( all_strands.End());
          strand_itr != strand_itr_end; ++strand_itr
        )
        {
          // increment sse_ctr
          ++strand_ctr;

          BCL_MessageStd( "Iteration #" + util::Format()( strand_ctr));

          // make a copy of the strand to be placed
          util::ShPtr< assemble::SSE> this_strand( ( *strand_itr)->HardCopy());
          this_strand->SetToIdealConformationAtOrigin();

          // call the Place and store the result
          storage::Pair< math::TransformationMatrix3D, bool> this_transformation
          (
            placement_always_flip.Place( *this_strand, iterative_model)
          );

          // apply the transformation and add it to model
          this_strand->Transform( this_transformation.First());
          iterative_model.Insert( this_strand);

          // collect the Sheets back from the model
          util::ShPtrVector< assemble::Domain> these_sheets( assemble::CollectorSheet().Collect( iterative_model));

          // make sure there is only one sheet of correct size
          BCL_MessageStd
          (
            " Number SSEs# " + util::Format()( iterative_model.GetNumberSSEs()) +
            " #sheets: " + util::Format()( these_sheets.GetSize()) + "\n" +
            these_sheets.FirstElement()->GetTopology()->GetOrderedIdentification()
          );

          // write the model
          Proteins::WriteModelToPDB
          (
            iterative_model,
            AddExampleOutputPathToFilename
            (
              placement, "placement_strand_next_to_sheet_itr_" + util::Format()( strand_ctr) + ".pdb"
            )
          );
        }

        // collect the Sheets back from the model
        util::ShPtrVector< assemble::Domain> these_sheets( assemble::CollectorSheet().Collect( iterative_model));

        BCL_Example_Check
        (
          these_sheets.GetSize() == 1 && these_sheets( 0)->GetNumberSSEs() == 5,
          "The end product should have been a single sheet with 5 strands but it is not, check the placement protocol!"
        );

      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldPlacementStrandNextToSheet

  const ExampleClass::EnumType ExampleFoldPlacementStrandNextToSheet::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldPlacementStrandNextToSheet())
  );

} // namespace bcl
