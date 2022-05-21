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
#include "math/bcl_math_roc_curve.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math_contingency_matrix.h"
#include "math/bcl_math_contingency_matrix_measures.h"
#include "math/bcl_math_piecewise_function.h"
#include "math/bcl_math_polynomial.h"
#include "math/bcl_math_quadratic_function.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_roc_curve.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathROCCurve :
    public ExampleInterface
  {
  public:

    ExampleMathROCCurve *Clone() const
    {
      return new ExampleMathROCCurve( *this);
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

      BCL_MessageStd( "Creating data set for classified and unclassified data");

      // initialize threshold
      const double threshold( 0.55);

      // create storage list of values and their actual unclassified expected values
      storage::List< storage::Pair< double, double> > values_unclassified;
      values_unclassified.PushBack( storage::Pair< double, double>( 0.00, 0.0));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.05, 0.1));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.06, 0.7));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.12, 0.6));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.15, 0.5));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.20, 0.4));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.26, 0.1));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.35, 0.1));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.38, 0.6));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.42, 0.3));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.47, 0.4));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.55, 0.8));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.63, 0.7));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.71, 0.3));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.75, 0.6));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.80, 0.7));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.82, 0.5));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.85, 0.8));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.90, 0.8));
      values_unclassified.PushBack( storage::Pair< double, double>( 0.99, 0.9));

      // create storage list of values and their actual classifications
      storage::List< storage::Pair< double, bool> > values_classified;
      values_classified.PushBack( storage::Pair< double, bool>( 0.00, false));
      values_classified.PushBack( storage::Pair< double, bool>( 0.05, false));
      values_classified.PushBack( storage::Pair< double, bool>( 0.06, true));
      values_classified.PushBack( storage::Pair< double, bool>( 0.12, true));
      values_classified.PushBack( storage::Pair< double, bool>( 0.15, false));
      values_classified.PushBack( storage::Pair< double, bool>( 0.20, false));
      values_classified.PushBack( storage::Pair< double, bool>( 0.26, false));
      values_classified.PushBack( storage::Pair< double, bool>( 0.35, false));
      values_classified.PushBack( storage::Pair< double, bool>( 0.38, true));
      values_classified.PushBack( storage::Pair< double, bool>( 0.42, false));
      values_classified.PushBack( storage::Pair< double, bool>( 0.47, false));
      values_classified.PushBack( storage::Pair< double, bool>( 0.55, true));
      values_classified.PushBack( storage::Pair< double, bool>( 0.63, true));
      values_classified.PushBack( storage::Pair< double, bool>( 0.71, false));
      values_classified.PushBack( storage::Pair< double, bool>( 0.75, true));
      values_classified.PushBack( storage::Pair< double, bool>( 0.80, true));
      values_classified.PushBack( storage::Pair< double, bool>( 0.82, false));
      values_classified.PushBack( storage::Pair< double, bool>( 0.85, true));
      values_classified.PushBack( storage::Pair< double, bool>( 0.90, true));
      values_classified.PushBack( storage::Pair< double, bool>( 0.99, true));

      // initialize expected integral
      const double expected_integral( 0.709695);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a ROC curve from classified data
      BCL_MessageStd( "Creating a ROC curve from classified data");
      // note: the list needs to be sorted previously - from large to small
      values_classified.Sort
      (
        std::binary_negate< storage::PairBinaryPredicateFirst< double, bool> >
        (
          storage::PairBinaryPredicateFirst< double, bool>
          (
            util::BinaryFunctionSTLWrapper< std::less< double> >()
          )
        )
      );
      math::ROCCurve roc_from_classified( values_classified);

      // create a ROC curve from unclassified data
      BCL_MessageStd( "Creating a ROC curve from unclassified data with a threshold of 0.55");
      // note: the list needs to be sorted previously - from large to small
      values_unclassified.Sort
      (
        std::binary_negate< storage::PairBinaryPredicateFirst< double, double> >
        (
          storage::PairBinaryPredicateFirst< double, double>
          (
            util::BinaryFunctionSTLWrapper< std::less< double> >()
          )
        )
      );
      math::ROCCurve roc_from_unclassified( values_unclassified, threshold);

    ////////////////
    // operations //
    ////////////////

      // create an empty ROC curves and then initialize them
      BCL_MessageStd( "Creating empty ROC curve and then initializing with classified data");
      math::ROCCurve roc_from_classified_b;
      roc_from_classified_b.Initialize( values_classified);

      BCL_MessageStd
      (
        "calculating integrals for roc curves from classified data and unclassified data"
      );

      // calculate integral
      const double integral_classified( roc_from_classified.Integral());

      // output integral
      BCL_MessageStd( "classified integral: " + util::Format()( integral_classified));

      // assert that classified integral is correct
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_integral, integral_classified, 0.001),
        "The integral from classified values should be equal to " + util::Format()( expected_integral)
        + " but is " + util::Format()( integral_classified)
      );

      // create an empty ROC curves and then initialize them
      BCL_MessageStd( "Creating empty ROC curve and then initializing with unclassified data");
      math::ROCCurve roc_from_unclassified_b( values_unclassified, threshold);

      // calculate integral
      const double integral_unclassified( roc_from_unclassified.Integral());

      // output integral
      BCL_MessageStd( "unclassified integral: " + util::Format()( integral_unclassified));

      // assert that unclassified integral is correct
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_integral, integral_unclassified, 0.001),
        "The integral from unclassified values should be equal to " + util::Format()( expected_integral)
        + " but is " + util::Format()( integral_classified)
      );

      // assert that fraction classified integral is correct
      BCL_ExampleCheckWithinTolerance( roc_from_classified_b.Integral( 0.5), 0.438596, 0.001);

      // check whether integral of roc curve can be reproduced by calling the generalized Integral method
      // using specific contigency function calls

      // check output integral of TPR vs FPR compared to ROC curve
      BCL_ExampleCheckWithinTolerance
      (
        roc_from_classified_b.Integral
        (
          double( 0.5),
          &math::ContingencyMatrix::GetFalsePositiveRate, // x-coord
          &math::ContingencyMatrix::GetTruePositiveRate   // y-coord
        ),
        0.438596,
        0.0001
      );

      // calculate weighted integral
      const double integral_unclassified_weighted( roc_from_unclassified.Integral( math::LinearFunction( -1, 1)));

      // output fraction weighted integral
      BCL_MessageStd
      (
        "the entire integral of classified weighted with f(x)=1-x : " + util::Format()( integral_unclassified_weighted)
      );

      // weighting function
      math::Polynomial made_polynomial( math::Polynomial::MakeFromCoefficients( linal::MakeVector< double>( 1, -2, 1)));

      // calculate weighted integral
      const double integral_unclassified_weighted_sqared( roc_from_unclassified.Integral( made_polynomial));

      // output fraction weighted integral
      BCL_MessageStd
      (
        "the entire integral of classified weighted with f(x)=(1-x)^2 : "
        + util::Format()( integral_unclassified_weighted_sqared)
      );

      // calculate integral of roc curve plotting FractionPredictedPositives (FPP) vs Precision (PPV)
      const double integral_ppv_vs_fpp
      (
        roc_from_classified_b.Integral
        (
          &math::ContingencyMatrix::GetFractionPredictedPositives, // x-coord
          &math::ContingencyMatrix::GetPrecision,                  // y-coord
          false
        )
      );

      // output
      BCL_MessageStd
      (
        "integral of roc curve plotting FractionPredictedPositives (FPP) vs Precision (PPV): "
        + util::Format()( integral_ppv_vs_fpp)
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( 0.72376, integral_ppv_vs_fpp, 0.001),
        "The value of integral calling specific x- and y- functions and the standard integral for roc curves should "
        "have equal values of 0.72376, but integral_generic is "
        + util::Format()( integral_ppv_vs_fpp)
      );

      // calculate integral of roc curve plotting FractionPredictedPositives (FPP) vs Precision (PPV)
      const double integral_ppv_vs_fpp_40pct
      (
        roc_from_classified_b.Integral
        (
          &math::ContingencyMatrix::GetFractionPredictedPositives, // x-coord
          &math::ContingencyMatrix::GetPrecision,                  // y-coord
          math::Range< double>( 0.0, 0.4)
        )
      );

      // output
      BCL_MessageStd
      (
        "integral of roc curve plotting FractionPredictedPositives (FPP) vs Precision (PPV) for 40pct of integral: "
        + util::Format()( integral_ppv_vs_fpp_40pct)
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( 0.881626, integral_ppv_vs_fpp_40pct, 0.001),
        "The value of integral calling specific x- and y- functions for 40pct of integral for roc curves should "
        "have equal values of 0.881626, but integral_generic is "
        + util::Format()( integral_ppv_vs_fpp_40pct)
      );

    /////////////////
    // data access //
    /////////////////

      // get the SortedCounts and output it
      BCL_MessageStd( "Outputting the sorted counts");
      BCL_MessageStd( util::Format()( roc_from_unclassified.GetSortedCounts()));

      // get the true negatives true positives and total number of results and output them
      BCL_MessageStd
      (
        "This roc curve consists of : " +
         util::Format()( roc_from_unclassified.GetNumberActualNegatives()) + "false positives " +
         util::Format()( roc_from_unclassified.GetNumberActualPositives()) + "true positives making " +
         util::Format()( roc_from_unclassified.GetNumberResults()) + "results in total"
      );

      // assert sum of true - true + = results
      BCL_Example_Check
      (
        roc_from_unclassified.GetNumberActualNegatives() +
        roc_from_unclassified.GetNumberActualPositives() ==  roc_from_unclassified.GetNumberResults(),
        " The sum of false positives and true positives does not match total number of results!!"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // initialize fstreams
      io::IFStream read;
      io::OFStream write;

      // outputting the ROC curve to a file
      BCL_MessageStd( "Writing the ROC curve to file and read it back in");
      WriteBCLObject( roc_from_unclassified);
      math::ROCCurve roc_curve_from_file;
      ReadBCLObject( roc_curve_from_file);

      // check the read ROC curve values
      BCL_MessageStd( "Checking the values for read roc curve\' values");
      const double integral_from_read_roc( roc_curve_from_file.Integral());

      BCL_MessageStd( "integral of read roc curve: " + util::Format()( integral_from_read_roc));

      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_integral, integral_from_read_roc, 0.001),
        "The integral from read roc curve should be " + util::Format()( expected_integral)
        + "but is " + util::Format()( integral_from_read_roc)
      );

      // outputting the ROC curve to a file for plotting
      BCL_MessageStd( "Writing the ROC curve to plot file test.roc_plot");
      std::string out_filename( AddExampleOutputPathToFilename( roc_from_unclassified, "test.roc_plot"));
      BCL_ExampleMustOpenOutputFile( write, out_filename);
      roc_from_unclassified.WritePlottingTable( write);
      io::File::CloseClearFStream( write);

      // read a ROC curve from the file
      BCL_MessageStd( "Read the ROC curve from plot file test.roc_plot");
      BCL_ExampleMustOpenInputFile( read, out_filename);
      math::ROCCurve roc_curve_from_plot_file;
      roc_curve_from_plot_file.ReadPlottingTable( read);
      io::File::CloseClearFStream( read);

      // check the read ROC curve values
      BCL_MessageStd( "Comparing the read roc plot files\' values");
      const double integral_from_read_roc_plot( roc_curve_from_plot_file.Integral());

      BCL_MessageStd
      (
        "integral of read roc plot file: " + util::Format()( integral_from_read_roc_plot)
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( double( 0.703324), integral_from_read_roc_plot, 0.001),
        "The integral from read roc plot file should be" + util::Format()( double( 0.703324))
        + " but is " + util::Format()( integral_from_read_roc_plot)
      );

      // outputting the ROC curve to a file for plotting after thinning
      BCL_MessageStd
      (
        "Writing a thinned (periodicity of 2) ROC curve to plot file test.roc_plot_thinned"
      );
      out_filename = AddExampleOutputPathToFilename( roc_from_unclassified, "test.roc_plot_thinned");
      BCL_ExampleMustOpenOutputFile( write, out_filename);
      roc_from_unclassified.GetThinnedRocCurvePeriodicity( 2).WritePlottingTable( write);
      io::File::CloseClearFStream( write);

      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( roc_from_unclassified, "test.roc_plot_thinned_formated"));
      roc_from_unclassified.GetThinnedRocCurvePeriodicity( 2).WriteRatePlottingTable( write, util::Format());
      io::File::CloseClearFStream( write);

      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( roc_from_unclassified, "test.roc_plot_ppv_vs_fpp_formated"));
      roc_from_unclassified.WriteRatePlottingTableGeneric
      (
        write,
        &math::ContingencyMatrix::GetFractionPredictedPositives, // x-coord
        &math::ContingencyMatrix::GetPrecision,                  // y-coord
        math::Range< double>( 0.0, 1.0),                         // x-coord cutoff
        util::Format(),
        std::string( "FractionPositivePrediced vs Precision")
      );
      io::File::CloseClearFStream( write);

      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( roc_from_unclassified, "test.roc_plot_ppv_vs_fpp_formated_40pct"));
      roc_from_unclassified.WriteRatePlottingTableGeneric
      (
        write,
        &math::ContingencyMatrix::GetFractionPredictedPositives, // x-coord
        &math::ContingencyMatrix::GetPrecision,                  // y-coord
        math::Range< double>( 0.0, 0.4),                         // x-coord cutoff range
        util::Format(),
        std::string( "FractionPositivePrediced vs Precision")
      );
      io::File::CloseClearFStream( write);

      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( roc_from_unclassified, "test.roc_plot_ppv_vs_fpp_formated_20_to_40pct"));
      roc_from_unclassified.WriteRatePlottingTableGeneric
      (
        write,
        &math::ContingencyMatrix::GetFractionPredictedPositives, // x-coord
        &math::ContingencyMatrix::GetPrecision,                  // y-coord
        math::Range< double>( 0.2, 0.4),                         // x-coord cutoff range
        util::Format(),
        std::string( "FractionPositivePrediced vs Precision")
      );
      io::File::CloseClearFStream( write);

      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( roc_from_unclassified, "test.roc_plot_ppv_vs_fpp_formated_20_to_40pct"));
      roc_from_unclassified.WriteRatePlottingTableGeneric
      (
        write,
        &math::ContingencyMatrix::GetFractionPredictedPositives, // x-coord
        &math::ContingencyMatrix::GetPrecision,                  // y-coord
        math::Range< double>( 0.2, 0.4),                         // x-coord cutoff range
        util::Format(),
        std::string( "FractionPositivePrediced vs Precision")
      );
      io::File::CloseClearFStream( write);
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( roc_from_unclassified, "test.roc_plot_maximized"));
      write << roc_from_unclassified.GetLocalPPVCurve();
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathROCCurve

  const ExampleClass::EnumType ExampleMathROCCurve::s_Instance
  (
    GetExamples().AddEnum( ExampleMathROCCurve())
  );

} // namespace bcl
