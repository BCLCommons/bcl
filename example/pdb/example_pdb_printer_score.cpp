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
#include "pdb/bcl_pdb_printer_score.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "fold/bcl_fold_default_scores.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_printer_score.cpp
  //!
  //! @author weinerbe
  //! @date Nov 18, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbPrinterScore :
    public ExampleInterface
  {
  public:

    ExamplePdbPrinterScore *Clone() const
    {
      return new ExamplePdbPrinterScore( *this);
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

      // get the default scores
      fold::DefaultScores::GetInstance().InitializeScores();
      fold::ScoreWeightSet score_weight_set;
      fold::DefaultScores::GetInstance().ModifyScoreWeightSet( score_weight_set);
      const util::ShPtr< score::ProteinModelScoreSum> sp_score( score_weight_set.ConstructScoreSum());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      const pdb::PrinterScore def_construct;

      // comment this out as currently the completeness_estimate is always printed when no quality measure is defined
      // so IsEmpty() is not true anymore.
//      util::GetLogger() << "DEBUG lines=" << def_construct( protein_model) << "\n";
//      BCL_ExampleIndirectCheck( def_construct( protein_model).IsEmpty(), true, "default construction");

      // test constructor from score
      const pdb::PrinterScore score_construct
      (
        sp_score,
        storage::Set< quality::Measure>( quality::GetMeasures().e_RMSD)
      );

    ///////////////
    // operators //
    ///////////////

      // test () operator
      const util::ShPtrList< pdb::Line> lines( score_construct( protein_model));
      BCL_ExampleCheck( lines.IsEmpty(), false);

      // write out the lines
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( def_construct, "score_printer.pdb"));
      pdb::Handler handler;
      handler.AppendLines( lines);
      handler.WriteLines( write);
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbPrinterScore

  const ExampleClass::EnumType ExamplePdbPrinterScore::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbPrinterScore())
  );

} // namespace bcl
