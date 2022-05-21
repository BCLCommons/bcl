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
#include "math/bcl_math_histogram.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_histogram.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathHistogram :
    public ExampleInterface
  {
  public:

    ExampleMathHistogram *Clone() const
    { return new ExampleMathHistogram( *this);}

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

      // this part is for the simple histogram
      const double values[ 20] = { 0, .5, .99, 7.0, 7.49999, 7.50000, 7.500001, 7.9999, 8, 9,
                            10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
      const storage::Vector< double> test( 20, values);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      math::Histogram histogram_default;

      BCL_Example_Check
      (
        ( histogram_default.GetBinSize()      == double( 0))
        && ( histogram_default.GetBoundaries().First() == double( 0))
        && ( histogram_default.GetBoundaries().Second() == double( 0))
        && ( histogram_default.GetNumberOfBins() == size_t( 0)),
        "constructing histogram from bin size, number bins and starting bin did not work"
      );

      // construct from starting bin, bin size and number of bins
      math::Histogram histogram1D_a( 0.7, 2.5, 10);
      math::Histogram histogram1D_b( 3.5, 0.7, 5);

      BCL_Example_Check
      (
        ( histogram1D_a.GetBinSize()      == double( 2.5))
        && ( histogram1D_a.GetBoundaries().First() == double( 0.7))
        && ( histogram1D_a.GetBoundaries().Second() == double( 0.7 + 10 * 2.5))
        && ( histogram1D_a.GetNumberOfBins() == size_t( 10)),
        "constructing histogram from bin size, number bins and starting bin did not work"
      );
      BCL_Example_Check
      (
        ( histogram1D_b.GetBinSize()      == double( 0.7))
        && ( histogram1D_b.GetBoundaries().First() == double( 3.5))
        && ( histogram1D_b.GetBoundaries().Second() == double( 3.5 + 5 * 0.7))
        && ( histogram1D_b.GetNumberOfBins() == size_t( 5)),
        "constructing histogram from bin size, number bins and starting bin did not work"
      );

      // clone
      util::ShPtr< math::Histogram> ptr( histogram1D_b.Clone());
      BCL_Example_Check
      (
        ( ptr->GetBinSize() == histogram1D_b.GetBinSize())
        && ( ptr->GetBoundaries() == histogram1D_b.GetBoundaries())
        && ( ptr->GetHistogram() == histogram1D_b.GetHistogram())
        && ( ptr->GetBoundariesCounts() == histogram1D_b.GetBoundariesCounts()),
        "clone did not copy the members correctly"
      );

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_MessageStd( "class name: " + ptr->GetClassIdentifier());
      BCL_Example_Check
      (
        GetStaticClassName< math::Histogram>() == "bcl::math::Histogram"
        && ptr->GetClassIdentifier() == GetStaticClassName< math::Histogram>(),
        "incorrect class name"
      );

      // access important parameters
      BCL_MessageStd( "binsize: "    + util::Format()( histogram1D_a.GetBinSize()));
      BCL_MessageStd( "boundaries: " + util::Format()( histogram1D_a.GetBoundaries()));
      BCL_MessageStd( "nr bins: "    + util::Format()( histogram1D_a.GetNumberOfBins()));

      // access the data
      BCL_MessageStd
      (
        "boundary counts: " + util::Format()( histogram1D_a.GetBoundariesCounts())
      );
      BCL_Example_Check
      (
        ( histogram1D_a.GetBoundariesCounts().First() == double( 0.0))
        && ( histogram1D_a.GetBoundariesCounts().Second() == double( 0.0)),
        "boundary counts should be 0.0, but are: " + util::Format()( histogram1D_a.GetBoundariesCounts())
      );

      // access the histogram counts
      BCL_Example_Check
      (
        ( histogram1D_a.GetHistogram() == linal::Vector< double>( 10, 0.0)),
        "histogram counts should be 0.0, but are: " + util::Format()( histogram1D_a.GetHistogram())
      );

      const double binning[] = { 1.95, 4.45, 6.95, 9.45, 11.95, 14.45, 16.95, 19.45, 21.95, 24.45};
      // the binning as vector - each number represents the middle of a bin
      BCL_Example_Check
      (
        ( histogram1D_a.GetBinning() == linal::Vector< double>( 10, binning)),
        "histogram binning should be " + util::Format()( linal::Vector< double>( 10, binning)) + " but is: " +
        util::Format()( histogram1D_a.GetBinning())
      );

      // the sum of all counts
      BCL_MessageStd
      (
        "sum counts: " + util::Format()( histogram1D_a.GetSumOfAllCounts())
      );
      BCL_Example_Check
      (
        ( histogram1D_a.GetSumOfAllCounts() == double( 0)),
        "histogram sum of all counts should be 0 but is: " + util::Format()( histogram1D_a.GetSumOfAllCounts())
      );

      // get the counts between two values
      BCL_MessageStd
      (
        "counts between 0 and 10: " +
        util::Format()( histogram1D_a.GetCountsInBetween( double( 0), double( 10)))
      );
      BCL_Example_Check
      (
        ( histogram1D_a.GetSumOfAllCounts() == double( 0)),
        "histogram counts between 0 and 10 should be 0 but are: " +
        util::Format()( histogram1D_a.GetCountsInBetween( double( 0), double( 10)))
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // calculate the histogram by passing in a vector of double values
      histogram1D_a.CalculateHistogram( test);

      // check the total count
      BCL_ExampleCheck( histogram1D_a.GetSumOfAllCounts(), double( 20));

      // check counts between 0 and 10
      BCL_ExampleCheck( histogram1D_a.GetCountsInBetween( 0.0, 10.0), double( 9));

      // check count for bin 2
      BCL_ExampleCheck( histogram1D_a.GetHistogram()( 2), double( 6));

      // PushBack an individual value - count in bin 2 should increase
      histogram1D_a.PushBack( double( 6.95));
      BCL_ExampleIndirectCheck( histogram1D_a.GetHistogram()( 2), double( 7), "histogram1D_a.PushBack( 6.95)");

      // PushBack an individual value of 0.5 - count in bin 2 should increase by 0.5
      histogram1D_a.PushBack( double( 6.95), double( 0.5));
      BCL_ExampleIndirectCheck( histogram1D_a.GetHistogram()( 2), double( 7.5), "histogram1D_a.PushBack( 6.95, 0.5)");

      // reset
      BCL_ExampleCheck( histogram1D_a.IsEmpty(), false);
      histogram1D_a.Reset();
      BCL_ExampleIndirectCheck( histogram1D_a.IsEmpty(), true, "histogram1D_a.Reset()");

      // last information containing bin
      BCL_ExampleCheck( math::Histogram().GetIndexOfLastInformationContainingBin(), 0);
      histogram1D_a.PushBack( double( 6.95), double( 0.5));
      BCL_ExampleCheck( histogram1D_a.GetIndexOfLastInformationContainingBin(), 2);

      // remove bins after index
      BCL_Example_Check
      (
        histogram1D_a.GetBoundariesCounts()( 1) == double( 0) && histogram1D_a.GetNumberOfBins() == 10,
        "boundary count should be 0 and number bins 10 but are: " +
        util::Format()( histogram1D_a.GetBoundariesCounts()( 1)) + " " +
        util::Format()( histogram1D_a.GetNumberOfBins())
      );
      histogram1D_a.RemoveBinsAfterIndex( 1);
      BCL_ExampleIndirectCheck( histogram1D_a.GetBoundariesCounts()( 1), 0.5, "histogram1D_a.RemoveBinsAfterIndex( 1)");
      BCL_ExampleIndirectCheck( histogram1D_a.GetNumberOfBins(), 2, "histogram1D_a.RemoveBinsAfterIndex( 1)");

      // combine two histograms
      math::Histogram combined_hist( histogram1D_a);
      const bool combine_aa_result( combined_hist.Combine( histogram1D_a));
      BCL_Example_Check
      (
        combine_aa_result
        && ( double( 2) * histogram1D_a.GetBoundariesCounts()( 0)) == combined_hist.GetBoundariesCounts()( 0)
        && ( double( 2) * histogram1D_a.GetBoundariesCounts()( 1)) == combined_hist.GetBoundariesCounts()( 1)
        && ( double( 2) * histogram1D_a.GetHistogram()) == combined_hist.GetHistogram(),
        "combined histogram should be twice the histogram1D_a"
      );

      // check combining with non identical histogram
      BCL_Example_Check
      (
        !histogram1D_a.Combine( histogram1D_b),
        "should not be able to combine histograms of different property"
      );

      // add pseudo count
      combined_hist.AddPseudoCount( double( 1));
      BCL_Example_Check
      (
        ( double( 2) * histogram1D_a.GetBoundariesCounts()( 0)) + double( 1) == combined_hist.GetBoundariesCounts()( 0)
        && ( double( 2) * histogram1D_a.GetBoundariesCounts()( 1)) + double( 1) == combined_hist.GetBoundariesCounts()( 1)
        && ( double( 2) * histogram1D_a.GetHistogram()) + double( 1) == combined_hist.GetHistogram(),
        "adding pseudo count did not work"
      );

      // prepare
      histogram1D_b = math::Histogram( 0.7, 2.5, 10);
      histogram1D_b.CalculateHistogram( test);

      // normalize
      histogram1D_b.Normalize();
      const double normalize_sum
      (
        histogram1D_b.GetBoundariesCounts().First() +
        histogram1D_b.GetHistogram().Sum() +
        histogram1D_b.GetBoundariesCounts().Second()
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( double( 1), normalize_sum),
        "normalization should result sum 1 but is: " + util::Format()( normalize_sum)
      );

      // calculate mean and sd
      const double approx_mean_b( math::Statistics::Mean( test.Begin(), test.End()));
      const double approx_sd_b( math::Statistics::StandardDeviation( test.Begin(), test.End()));
      BCL_MessageStd( "mean and sd from unbinned data: " + util::Format()( approx_mean_b) + " " + util::Format()( approx_sd_b));

      const double expected_mean_b( 10.9778);
      const double expected_sd_b( 4.65019);
      const double calculated_mean_b( histogram1D_b.CalculateMean());
      const double calculated_sd_b( histogram1D_b.CalculateSD());
      BCL_MessageStd( "mean and sd calculated from histogram: " + util::Format()( calculated_mean_b) + " " + util::Format()( calculated_sd_b));

      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_mean_b, calculated_mean_b),
        "calc and expected mean do not agree " + util::Format()( calculated_mean_b) + " != " + util::Format()( expected_mean_b)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_sd_b, calculated_sd_b),
        "calc and expected sd do not agree " + util::Format()( calculated_sd_b) + " != " + util::Format()( expected_sd_b)
      );

      // ExtendBoundaries
      {
        {
          math::Histogram short_histogram( 15.0, 1, 10);
          BCL_ExampleCheck( true, short_histogram.ExtendBoundaries( 2, 0, 13, 3, 0, 28));
          BCL_ExampleCheck( short_histogram.GetNumberOfBins(), 15);
          BCL_ExampleCheck( short_histogram.GetBoundaries().First(), 13);
          BCL_ExampleCheck( short_histogram.GetBoundaries().Second(), 28);
        }
        {
          math::Histogram short_histogram( 15.0, 1, 10);
          BCL_ExampleCheck( true, short_histogram.ExtendBoundaries( 0, 0, 15, 2, 0, 27));
          BCL_ExampleCheck( short_histogram.GetNumberOfBins(), 12);
          BCL_ExampleCheck( short_histogram.GetBoundaries().First(), 15);
          BCL_ExampleCheck( short_histogram.GetBoundaries().Second(), 27);
        }
        {
          math::Histogram short_histogram( 15.0, 1, 10);
          BCL_ExampleCheck( true, short_histogram.ExtendBoundaries( 4, 0, 11, 0, 0, 25));
          BCL_ExampleCheck( short_histogram.GetNumberOfBins(), 14);
          BCL_ExampleCheck( short_histogram.GetBoundaries().First(), 11);
          BCL_ExampleCheck( short_histogram.GetBoundaries().Second(), 25);
        }
        {
          math::Histogram short_histogram( 15.0, 1, 10);
          BCL_ExampleCheck( true, short_histogram.ExtendBoundaries( 4, 0, 12, 0, 0, 25));
          BCL_ExampleCheck( short_histogram.GetNumberOfBins(), 13);
          BCL_ExampleCheck( short_histogram.GetBoundaries().First(), 12);
          BCL_ExampleCheck( short_histogram.GetBoundaries().Second(), 25);
        }
        {
          math::Histogram short_histogram( 15.0, 1, 10);
          BCL_ExampleCheck( true, short_histogram.ExtendBoundaries( 0, 0, 15, 2, 0, 26));
          BCL_ExampleCheck( short_histogram.GetNumberOfBins(), 11);
          BCL_ExampleCheck( short_histogram.GetBoundaries().First(), 15);
          BCL_ExampleCheck( short_histogram.GetBoundaries().Second(), 26);
        }
        {
          math::Histogram short_histogram( 15.0, 1, 10);
          BCL_ExampleCheck( true, short_histogram.ExtendBoundaries( 2, 0, 13.3, 3, 0, 27.6));
          BCL_ExampleCheck( short_histogram.GetNumberOfBins(), 13);
          BCL_ExampleCheck( short_histogram.GetBoundaries().First(), 14);
          BCL_ExampleCheck( short_histogram.GetBoundaries().Second(), 27);
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

      // prepare
      histogram1D_b = math::Histogram( 0.7, 2.5, 10);
      histogram1D_b.CalculateHistogram( test);

      // write to file
      WriteBCLObject( histogram1D_b);
      // read
      ReadBCLObject( histogram_default);

      BCL_Example_Check
      (
        ( histogram_default.GetBinSize() == histogram1D_b.GetBinSize())
        && ( histogram_default.GetBoundaries() == histogram1D_b.GetBoundaries())
        && ( histogram_default.GetHistogram() == histogram1D_b.GetHistogram())
        && ( histogram_default.GetBoundariesCounts() == histogram1D_b.GetBoundariesCounts()),
        "read histogram is different from written histogram\n" + util::Format()( histogram1D_b) + "\n" +
        util::Format()( histogram_default)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleMathHistogram

  const ExampleClass::EnumType ExampleMathHistogram::s_Instance
  (
    GetExamples().AddEnum( ExampleMathHistogram())
  );

} // namespace bcl
