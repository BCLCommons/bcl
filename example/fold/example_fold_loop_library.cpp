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
#include "fold/bcl_fold_loop_library.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_loop_library.cpp
  //! @brief this example tests the implementation of fold::LoopLibrary
  //!
  //! @author fischea
  //! @date Dec 18, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleFoldLoopLibrary :
     public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief clone function
    //! @return pointer to a new ExampleFoldLoopLibrary
    ExampleFoldLoopLibrary *Clone() const
    {
      return new ExampleFoldLoopLibrary( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {
      // prepare variables to instantiate the loop template library
      const std::string lib_file_name( AddExampleInputPathToFilename( e_Fold, "loop_library.ls"));
      io::IFStream lib_file;
      BCL_ExampleMustOpenInputFile( lib_file, lib_file_name);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create test handler
      util::ShPtr< fold::LoopLibrary> sp_library( fold::LoopLibrary::CreateLoopLibrary( lib_file_name));

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( sp_library->GetClassIdentifier(), ( GetStaticClassName< fold::LoopLibrary>()));

    ////////////////
    // operations //
    ////////////////

      // prepare test files
      BCL_ExampleMustOpenInputFile( lib_file, lib_file_name);
      size_t num_loops;
      fold::LoopParameters loop_1;
      lib_file >> num_loops;
      lib_file >> loop_1;
      fold::LoopParameters loop_2;
      lib_file >> loop_2;
      io::File::CloseClearFStream( lib_file);

      // check if the correct template for the first loop is found
      const fold::LoopParameters loop_1_temp( *sp_library->FindTemplates( loop_1).FirstElement());
      BCL_ExampleCheck( loop_1_temp.GetRotation(), loop_1.GetRotation());
      BCL_ExampleCheck( loop_1_temp.GetTranslation(), loop_1.GetTranslation());
      BCL_ExampleCheck( loop_1_temp.GetSequenceDistance(), loop_1.GetSequenceDistance());

      // check if the correct template for the second loop is found
      const fold::LoopParameters loop_2_temp( *sp_library->FindTemplates( loop_2).FirstElement());
      BCL_ExampleCheck( loop_2_temp.GetRotation(), loop_2.GetRotation());
      BCL_ExampleCheck( loop_2_temp.GetTranslation(), loop_2.GetTranslation());
      BCL_ExampleCheck( loop_2_temp.GetSequenceDistance(), loop_2.GetSequenceDistance());

      // recombine templates to get a template of length 10
      const std::string rec_lib_file_name( AddExampleInputPathToFilename( e_Mc, "loop_recombination_test_lib"));
      util::ShPtr< fold::LoopLibrary> sp_rec_library( fold::LoopLibrary::CreateLoopLibrary( rec_lib_file_name));

      // identify the shorter templates as well as the combined template used as reference in the library
      storage::HashMap< std::string, util::ShPtrVector< fold::LoopParameters> > templates
      (
        sp_rec_library->GetTemplates()
      );
      util::ShPtr< fold::LoopParameters> sp_template_n;
      util::ShPtr< fold::LoopParameters> sp_template_c;
      util::ShPtr< fold::LoopParameters> sp_template_combined;
      for( auto temp_it( templates.Begin()), temp_it_end( templates.End()); temp_it != temp_it_end; ++temp_it)
      {
        const util::ShPtr< fold::LoopParameters> sp_template_tmp( temp_it->second( 0));
        if( sp_template_tmp->GetSequenceDistance() == 2)
        {
          if( !sp_template_c.IsDefined())
          {
            sp_template_c = sp_template_tmp;
          }
          {
            sp_template_n = sp_template_tmp;
          }
        }
        else
        {
          sp_template_combined = sp_template_tmp;
        }
      }

      // check if templates can be sampled
      const util::ShPtrVector< fold::LoopParameters> new_templates( sp_rec_library->GenerateTemplates( 5, 1));
      BCL_Example_Check
      (
        new_templates.GetSize() == 1 && new_templates( 0)->GetSequenceDistance() == 5,
        "generate templates from existing templates"
      );

      // check if templates can be sampled
      storage::HashMap< size_t, size_t> counts;
      counts.Insert( storage::Pair< size_t, size_t>( 5, 1));
      const size_t old_number( sp_rec_library->FindTemplates( 5).GetSize());
      sp_rec_library->RecombineTemplates( counts);
      const size_t new_number( sp_rec_library->FindTemplates( 5).GetSize());
      BCL_ExampleCheck( new_number, old_number + 1);

      // check template requirement estimation
      const std::string model_filename( AddExampleInputPathToFilename( e_Fold, "1x91_no_loops.pdb"));
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      const util::ShPtr< assemble::ProteinModel> sp_model
      (
        util::CloneToShPtr( factory.ProteinModelFromPDBFilename( model_filename))
      );
      util::ShPtrVector< assemble::ProteinModel> models;
      models.PushBack( sp_model);
      const util::ShPtr< storage::HashMap< size_t, size_t> > sp_estimate
      (
        fold::LoopLibrary::EstimateTemplateRequirement( models)
      );

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

  }; // class ExampleFoldLoopLibrary

  //! single instance of this class
  const ExampleClass::EnumType ExampleFoldLoopLibrary::s_Instance
  (
     GetExamples().AddEnum( ExampleFoldLoopLibrary())
  );

} // namespace bcl
