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
#include "fold/bcl_fold_mutate_protein_model_compress.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_compress.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelCompress :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelCompress *Clone() const
    {
      return new ExampleFoldMutateProteinModelCompress( *this);
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

      //build models from pdb
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> sse_min_sizes;
      sse_min_sizes[ biol::GetSSTypes().HELIX] = 9;
      sse_min_sizes[ biol::GetSSTypes().STRAND] = 5;
      sse_min_sizes[ biol::GetSSTypes().COIL] = 999;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, sse_min_sizes)
      );

      // initialize compression factors to be used
      const math::Range< double> compression_range( 0.95, 1.02);
      const double compression_factor( 0.95);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      fold::MutateProteinModelCompress mutate_default;

      // test construction from a compression factor
      BCL_MessageStd( "test construct from a compression factor");
      fold::MutateProteinModelCompress mutate_factor( compression_factor);

      // test construction from a compression range
      BCL_MessageStd( "test construct from a compression range");
      fold::MutateProteinModelCompress mutate_range( compression_range);

      // test copy constructor
      fold::MutateProteinModelCompress mutate_copy( mutate_range);
      BCL_MessageStd( "test copy constructor");
      BCL_Example_Check
      (
        mutate_copy.GetCompressionRange() == mutate_range.GetCompressionRange(),
        "copy should be " + util::Format()( mutate_range) + " but instead is " + util::Format()( mutate_copy)
      );

      // test clone constructor
      util::ShPtr< fold::MutateProteinModelCompress> mutate_clone( mutate_copy.Clone());
      BCL_MessageStd( "test clone constructor");
      BCL_Example_Check
      (
        mutate_clone->GetCompressionRange() == mutate_copy.GetCompressionRange(),
        "clone  should be " + util::Format()( mutate_copy) + " but instead is " + util::Format()( *mutate_clone)
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      BCL_MessageStd( "test GetStatisClassName");
      const std::string correct_static_class_name( "bcl::fold::MutateProteinModelCompress");
      BCL_Example_Check
      (
        GetStaticClassName< fold::MutateProteinModelCompress>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< fold::MutateProteinModelCompress>() + " but should give " +
        correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_MessageStd( "test GetClassIdentifier");
      BCL_Example_Check
      (
        GetStaticClassName< fold::MutateProteinModelCompress>() == mutate_factor.GetClassIdentifier(),
        "GetClassIdentifier gives " + mutate_factor.GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // check GetCompressionfactorRange()
      BCL_MessageStd( "test GetCompressionfactorRange");
      BCL_Example_Check
      (
        mutate_range.GetCompressionRange() == compression_range,
        "GetCompressionfactorRange gives " + util::Format()( mutate_range.GetCompressionRange()) +
        " instead of " + util::Format()( compression_range)
      );

    ///////////////
    // operators //
    ///////////////

      // test operator =
      fold::MutateProteinModelCompress mutate_assign;
      mutate_assign = mutate_range;
      BCL_MessageStd( "test operator =");
      BCL_Example_Check
      (
        mutate_assign.GetCompressionRange() == mutate_range.GetCompressionRange(),
        "assigned should be " + util::Format()( mutate_range) + "\nbut instead is\n" + util::Format()( mutate_assign)
      );

      // initialize variables
      const util::SiPtrVector< const assemble::SSE> sses( protein_model.GetSSEs());
      linal::Vector3D center_of_mass( protein_model.GetCenterOfSSEs());
      linal::Vector< double> sse_center_distances( sses.GetSize());

      // initialize SSE ctr
      size_t sse_ctr( 0);
      // iterate over SSEs
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end; ++sse_itr, ++sse_ctr
      )
      {
        // update distance from this SSE center to center of mass
        sse_center_distances( sse_ctr) = linal::Distance( ( *sse_itr)->GetCenter(), center_of_mass);
      }

      BCL_MessageStd( "test operator() in an iteration");
      bool success_center_of_mass( true);
      bool success_sse_distances( true);

      // iterate 3 times
      for( size_t ctr( 0); ctr < 3; ++ctr)
      {
        BCL_MessageStd( "iteration #" + util::Format()( ctr));
        // apply the mutate and store the new model
        BCL_MessageStd( "applying mutation");
        protein_model = ( *mutate_factor( protein_model).GetArgument());

        // write the new model to file
        BCL_MessageStd( "writing model to file");
        const std::string
          output_filename( AddExampleOutputPathToFilename( mutate_default, "1ubi.pdb.compressed." + util::Format()( ctr)));
        Proteins::WriteModelToPDB( protein_model, output_filename);

        // calculate the new center of mass
        const linal::Vector3D new_center_of_mass( protein_model.GetCenterOfSSEs());

        // now check the new center of mass
        BCL_MessageStd( "checking new center of mass");
        BCL_MessageStd
        (
          "old center of mass :\n" + util::Format()( center_of_mass) + "\n" +
          "new center of mass :\n" + util::Format()( new_center_of_mass)
        );
        success_center_of_mass &= math::EqualWithinTolerance( center_of_mass, new_center_of_mass);

        // calculate SSE distances to center
        const util::SiPtrVector< const assemble::SSE> new_sses( protein_model.GetSSEs());
        linal::Vector< double> new_sse_center_distances( new_sses.GetSize());

        BCL_MessageStd( "Checking distances from SSE centers to new center of mass");
        // initialize SSE ctr
        size_t sse_ctr( 0);
        // iterate over SSEs
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( new_sses.Begin()),
            sse_itr_end( new_sses.End());
          sse_itr != sse_itr_end; ++sse_itr, ++sse_ctr
        )
        {
          // calculate distance from this SSE center to center of mass
          new_sse_center_distances( sse_ctr) = linal::Distance( ( *sse_itr)->GetCenter(), new_center_of_mass);

          // compare the value with the previous one
          BCL_MessageStd
          (
            ( *sse_itr)->GetIdentification() +
            " old: " + util::Format()( sse_center_distances( sse_ctr)) +
            " new: " + util::Format()( new_sse_center_distances( sse_ctr))
          );
          success_sse_distances &=
            math::EqualWithinTolerance( compression_factor * sse_center_distances( sse_ctr), new_sse_center_distances( sse_ctr));
        }

        // update values to recent ones for the next check
        center_of_mass = new_center_of_mass;
        sse_center_distances = new_sse_center_distances;
      }

      BCL_Example_Check
      (
        success_center_of_mass, "the center of mass has changed significantly"
      );

      BCL_Example_Check
      (
        success_sse_distances, "the distances differ significantly"
      );

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      WriteBCLObject( mutate_range);
      fold::MutateProteinModelCompress mutate_read;
      ReadBCLObject( mutate_read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_Example_Check
      (
        mutate_read.GetCompressionRange() == mutate_range.GetCompressionRange(),
        "the written and read MutateProteinModelCompress classes differ from each other" +
        util::Format()( mutate_read) + "\nvs\n" + util::Format()( mutate_range)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelCompress

  const ExampleClass::EnumType ExampleFoldMutateProteinModelCompress::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelCompress())
  );

} // namespace bcl
