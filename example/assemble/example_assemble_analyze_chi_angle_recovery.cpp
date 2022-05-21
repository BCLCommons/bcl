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
#include "assemble/bcl_assemble_analyze_chi_angle_recovery.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_rotamer.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_analyze_chi_angle_recovery.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Aug 27, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleAnalyzeChiAngleRecovery :
    public ExampleInterface
  {
  public:

    ExampleAssembleAnalyzeChiAngleRecovery *Clone() const
    {
      return new ExampleAssembleAnalyzeChiAngleRecovery( *this);
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
      const std::string native_chi_filename
      (
        AddExampleInputPathToFilename( e_Biology, "AnalyzeChiAngleRecovery.native_clone")
      );

      BCL_MessageDbg( "native_chi_filename is " + util::Format()( native_chi_filename));

      const std::string postfix( ".ChiAngleRecovery");

      const std::string label
      (
        "ChiAngleRecovery ( filename_postfix=" + postfix + ",collector_type=CollectorAAType( "
        "methanesulfonothioate ) ,"
        "reference_chi_filename=" + native_chi_filename + ",tolerance=10,"
        "angle_unit=degree )"
      );

      util::Implementation< assemble::AnalyzeProteinEnsembleInterface> analysis( label);

      const std::string ensemble_filename
      (
        AddExampleInputPathToFilename( e_Biology, "AnalyzeChiAnglePairDistribution_ensemble.ls")
      );

      const std::string example_prefix( AddExampleInputPathToFilename( e_Biology, ""));
      const assemble::ProteinEnsemble ensemble( ensemble_filename, 0, biol::GetAAClasses().e_AAComplete, example_prefix);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      assemble::AnalyzeChiAngleRecovery def_constr;

      // clone constructor
      {
        util::ShPtr< assemble::AnalyzeProteinEnsembleInterface> clone_constr( analysis->Clone());
        const std::string prefix( AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAngleRecovery_clone"));
        analysis->WriteAnalysisFile( prefix, ensemble);

        BCL_ExampleCheck
        (
          io::File::FilesMatch
          (
            AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAngleRecovery_clone_correct" + postfix),
            prefix + postfix
          ),
          true
        );
      }

    /////////////////
    // data access //
    /////////////////

      // GetOutFilePostfix
      BCL_ExampleCheck( analysis->GetOutFilePostfix(), postfix);

    ///////////////
    // operators //
    ///////////////

      // operator
      // check previous chi dependence - chi 4 is wrong but chi 5 is right so chi 1-3 should be correct
      {
        const std::string native_chi_filename
        (
          AddExampleInputPathToFilename( e_Biology, "AnalyzeChiAngleRecovery.native_b")
        );

        const std::string label
        (
          "ChiAngleRecovery ( filename_postfix=" + postfix + ",collector_type=CollectorAAType("
          "methanesulfonothioate) ,"
          "reference_chi_filename=" + native_chi_filename + ",tolerance=10,"
          "angle_unit=degree )"
        );

        util::Implementation< assemble::AnalyzeProteinEnsembleInterface> analysis( label);
        util::ShPtr< assemble::AnalyzeProteinEnsembleInterface> clone_constr( analysis->Clone());
        const std::string prefix( AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAngleRecovery_b"));
        analysis->WriteAnalysisFile( prefix, ensemble);

        BCL_ExampleCheck
        (
          io::File::FilesMatch
          (
            AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAngleRecovery_b_correct" + postfix),
            prefix + postfix
          ),
          true
        );
      }
      // check with nan - as soon as nan is seen following chis cant be correct
      {
        const std::string native_chi_filename
        (
          AddExampleInputPathToFilename( e_Biology, "AnalyzeChiAngleRecovery.native_c")
        );

        const std::string label
        (
          "ChiAngleRecovery ( filename_postfix=" + postfix + ",collector_type=CollectorAAType("
          "methanesulfonothioate ) ,"
          "reference_chi_filename=" + native_chi_filename + ",tolerance=10,"
          "angle_unit=degree )"
        );

        util::Implementation< assemble::AnalyzeProteinEnsembleInterface> analysis( label);
        util::ShPtr< assemble::AnalyzeProteinEnsembleInterface> clone_constr( analysis->Clone());
        const std::string prefix( AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAngleRecovery_c"));
        analysis->WriteAnalysisFile( prefix, ensemble);

        BCL_ExampleCheck
        (
          io::File::FilesMatch
          (
            AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAngleRecovery_c_correct" + postfix),
            prefix + postfix
          ),
          true
        );
      }
      // check with more than one native rotamer
      {
        const std::string native_chi_filename
        (
          AddExampleInputPathToFilename( e_Biology, "AnalyzeChiAngleRecovery.native_d")
        );

        const std::string label
        (
          "ChiAngleRecovery ( filename_postfix=" + postfix + ",collector_type=CollectorAAType("
          "methanesulfonothioate) ,"
          "reference_chi_filename=" + native_chi_filename + ",tolerance=15,"
          "angle_unit=degree )"
        );

        util::Implementation< assemble::AnalyzeProteinEnsembleInterface> analysis( label);
        util::ShPtr< assemble::AnalyzeProteinEnsembleInterface> clone_constr( analysis->Clone());
        const std::string prefix( AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAngleRecovery_d"));
        analysis->WriteAnalysisFile( prefix, ensemble);

        BCL_ExampleCheck
        (
          io::File::FilesMatch
          (
            AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAngleRecovery_d_correct" + postfix),
            prefix + postfix
          ),
          true
        );
      }
      // check with more than one native rotamer last rotamer does not have bcl::biol::Rotamer
      {
        const std::string native_chi_filename
        (
          AddExampleInputPathToFilename( e_Biology, "AnalyzeChiAngleRecovery.native_e")
        );

        const std::string label
        (
          "ChiAngleRecovery ( filename_postfix=" + postfix + ",collector_type=CollectorAAType("
          "methanesulfonothioate) ,"
          "reference_chi_filename=" + native_chi_filename + ",tolerance=15,"
          "angle_unit=degree )"
        );

        util::Implementation< assemble::AnalyzeProteinEnsembleInterface> analysis( label);
        util::ShPtr< assemble::AnalyzeProteinEnsembleInterface> clone_constr( analysis->Clone());
        const std::string prefix( AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAngleRecovery_e"));
        analysis->WriteAnalysisFile( prefix, ensemble);

        BCL_ExampleCheck
        (
          io::File::FilesMatch
          (
            AddExampleOutputPathToFilename( ensemble, "AnalyzeChiAngleRecovery_e_correct" + postfix),
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

      // read write unnecessary

    //////////////////////
    // helper functions //
    //////////////////////

      {
        assemble::ProteinModel model
        (
          pdb::Factory( biol::GetAAClasses().e_AAComplete).ProteinModelFromPDBFilename
          (
            AddExampleInputPathToFilename( e_Biology, "3mpq.pdb")
          )
        );

        assemble::CollectorAAType collector( biol::GetAATypes().R1A);
        const util::SiPtrList< const biol::AABase> aas( collector.Collect( model.GetAminoAcids()));
        BCL_Assert( aas.GetSize() == 1, "not exactly one aa collected");
        const biol::Rotamer rotamer( aas.FirstElement()->CalculateSideChainDihedralAngles());
        BCL_MessageDbg( "3mpn rotamer is " + rotamer.WriteSimple( math::Angle::e_Degree));
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleAnalyzeChiAngleRecovery

  const ExampleClass::EnumType ExampleAssembleAnalyzeChiAngleRecovery::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleAnalyzeChiAngleRecovery())
  );

} // namespace bcl
