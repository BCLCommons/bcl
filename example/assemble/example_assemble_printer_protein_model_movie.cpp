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
#include "assemble/bcl_assemble_printer_protein_model_movie.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "fold/bcl_fold_default_scores.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "mc/bcl_mc_movie_printer_pymol.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_printer_protein_model_movie.cpp
  //!
  //! @author weinerbe
  //! @date Nov 28, 2011
  //! @remarks status incomplete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssemblePrinterProteinModelMovie :
    public ExampleInterface
  {
  public:

    ExampleAssemblePrinterProteinModelMovie *Clone() const
    {
      return new ExampleAssemblePrinterProteinModelMovie( *this);
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
      // get a protein model
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );
      storage::Pair< util::ShPtr< assemble::ProteinModel>, double> model_pair( util::PtrHardCopy( &protein_model), 0.0);

      // initialize superimposition
      const quality::SuperimposeMeasure &superimpose( quality::GetSuperimposeMeasures().e_RMSD);

      // get the default scores
      fold::DefaultScores::GetInstance().InitializeScores();
      fold::ScoreWeightSet score_weight_set;
      fold::DefaultScores::GetInstance().ModifyScoreWeightSet( score_weight_set);
      const util::ShPtr< score::ProteinModelScoreSum> sp_score( score_weight_set.ConstructScoreSum());

      // initialize step status
      const storage::Set< opti::StepStatusEnum> step_status_set( opti::e_Improved);

      // initialize quality measures
      const storage::Set< quality::Measure> qualities( quality::GetMeasures().e_RMSD);

      // initialize movie pointer
      const std::string prefix( AddExampleOutputPathToFilename( assemble::GetNamespaceIdentifier(), "movie_printer_"));
      util::ShPtr< mc::MoviePrinterInterface> sp_movie( new mc::MoviePrinterPymol());
      sp_movie->Initialize
      (
        prefix,
        storage::Vector< std::string>::Create
        (
          GetStaticClassName< storage::Table< double> >(),
          math::SumFunctionMixin< score::ProteinModel>::GetValueTableVerticalColumnNames()( 2)
        ),
        sp_score->GetFunctionSchemes(),
        720,
        480,
        false
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::PrinterProteinModelMovie def_construct;
      BCL_ExampleCheck( def_construct.GetPrefix(), "");

      // construct from members
      assemble::PrinterProteinModelMovie movie_printer
      (
        prefix,
        sp_movie,
        sp_score,
        step_status_set,
        superimpose,
        qualities
      );

    /////////////////
    // data access //
    /////////////////

      // test GetPrefix
      BCL_ExampleCheck( movie_printer.GetPrefix(), prefix);

    ////////////////
    // operations //
    ////////////////

      // initialize
      movie_printer.Initialize( 13, 2);

      // test Print functions
      // TODO

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssemblePrinterProteinModelMovie

  const ExampleClass::EnumType ExampleAssemblePrinterProteinModelMovie::s_Instance
  (
    GetExamples().AddEnum( ExampleAssemblePrinterProteinModelMovie())
  );

} // namespace bcl
