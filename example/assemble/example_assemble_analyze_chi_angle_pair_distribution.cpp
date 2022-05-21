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
#include "assemble/bcl_assemble_analyze_chi_angle_pair_distribution.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_aa_classes.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_analyze_chi_angle_pair_distribution.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Aug 22, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleAnalyzeChiAnglePairDistribution :
    public ExampleInterface
  {
  public:

    ExampleAssembleAnalyzeChiAnglePairDistribution *Clone() const
    {
      return new ExampleAssembleAnalyzeChiAnglePairDistribution( *this);
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
      const std::string reference_pairs_filename
      (
        AddExampleInputPathToFilename( e_Biology, "AnalyzeChiAnglePairDistribution_reference_pairs.ls")
      );

      BCL_MessageDbg( "reference pairs filename is " + util::Format()( reference_pairs_filename));

      const std::string postfix( ".ChiAnglePairDistribution");

      const std::string label
      (
        "ChiAnglePairDistribution ( filename_postfix=" + postfix + ",histogram_binsize=10,"
        "collector_type=CollectorAAType( methanesulfonothioate ) ,chi_angle_a=e_One, chi_angle_b=e_Two, "
        "set_title_and_label=0,title=TitleTest, show_color_box=0,font_size=16,x_pixels=640,y_pixels=640,"
        "reference_chi_filename=" + reference_pairs_filename + " )"
      );

      util::Implementation< assemble::AnalyzeProteinEnsembleInterface> analysis( label);

      const std::string ensemble_filename
      (
        AddExampleInputPathToFilename( e_Biology, "AnalyzeChiAnglePairDistribution_ensemble.ls")
      );

      const std::string example_prefix( AddExampleInputPathToFilename( e_Biology, ""));
      const assemble::ProteinEnsemble ensemble( ensemble_filename, 0, biol::GetAAClasses().e_AAComplete, example_prefix);
      const std::string output( analysis->operator()( ensemble));
      BCL_MessageDbg( "output is \n" + output + "\n");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const assemble::AnalyzeChiAnglePairDistribution def_constr;

      // clone constructor
      util::ShPtr< assemble::AnalyzeProteinEnsembleInterface> clone_constr( analysis->Clone());
      {
        const std::string prefix( AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAnglePairDistribution_clone"));
        analysis->WriteAnalysisFile( prefix, ensemble);

        BCL_ExampleCheck
        (
          io::File::FilesMatch
          (
            AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAnglePairDistribution_correct" + postfix),
            prefix + postfix
          ),
          true
        );
      }

    /////////////////
    // data access //
    /////////////////

      // GetMinMax
      const storage::VectorND< 2, double> min_max( -180, 180);
      BCL_ExampleCheck( min_max, assemble::AnalyzeChiAnglePairDistribution::GetMinMax());

      // GetOutFilePostfix
      BCL_ExampleCheck( def_constr.GetOutFilePostfix(), postfix);

    ///////////////
    // operators //
    ///////////////

      {
        const std::string prefix( AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAnglePairDistribution"));
        analysis->WriteAnalysisFile( prefix, ensemble);

        BCL_ExampleCheck
        (
          io::File::FilesMatch
          (
            AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAnglePairDistribution_correct" + postfix),
            prefix + postfix
          ),
          true
        );
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test Write and read
      WriteBCLObject( *clone_constr);
      assemble::AnalyzeChiAnglePairDistribution read_analysis;
      ReadBCLObject( read_analysis);

      {
        const std::string prefix( AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAnglePairDistribution_read"));
        analysis->WriteAnalysisFile( prefix, ensemble);

        BCL_ExampleCheck
        (
          io::File::FilesMatch
          (
            AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAnglePairDistribution_correct" + postfix),
            prefix + postfix
          ),
          true
        );
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleAnalyzeChiAnglePairDistribution

  const ExampleClass::EnumType ExampleAssembleAnalyzeChiAnglePairDistribution::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleAnalyzeChiAnglePairDistribution())
  );

} // namespace bcl
