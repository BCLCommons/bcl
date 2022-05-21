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
#include "math/bcl_math_histogram_2d.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_histogram_2d.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathHistogram2D :
    public ExampleInterface
  {
  public:

    ExampleMathHistogram2D *Clone() const
    {
      return new ExampleMathHistogram2D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      size_t number_value_pairs( 10000);
      storage::VectorND< 2, double> minmax_x( 0, 3);
      storage::VectorND< 2, double> minmax_y( 0, 3);
      storage::VectorND< 2, double> binsize_xy( 0.5, 0.25);
      storage::VectorND< 2, size_t> number_of_bins_xy( 6, 12);
      storage::Vector< storage::VectorND< 2, double> > values_vector( number_value_pairs);

      //fill valuesvector with random values
      BCL_MessageStd
      (
        std::string( "filling a storage::Vector of Pairs of double with random valuepairs")
      );
      for( size_t i( 0); i < number_value_pairs; ++i)
      {
        values_vector( i)
          = storage::VectorND< 2, double>
            (
              random::GetGlobalRandom().Random< double>( minmax_x.First() - 2, minmax_x.Second() + 2),
              random::GetGlobalRandom().Random< double>( minmax_y.First() - 2, minmax_y.Second() + 2)
            );
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      math::Histogram2D histogram2d_default;

      //instantiate histogram2D
      BCL_MessageStd( std::string( "building and calculating histogram"));
      math::Histogram2D histogram
      (
        storage::VectorND< 2, double>( minmax_x.First(), minmax_y.First()), binsize_xy, number_of_bins_xy
      );
      BCL_Example_Check
      (
        ( histogram.GetBinSizeXY()            == binsize_xy)
        && ( histogram.GetBoundariesX().First()  == minmax_x.First())
        && ( histogram.GetBoundariesY().First()  == minmax_y.First())
        && ( histogram.GetBoundariesX().Second() == minmax_x.First() + number_of_bins_xy.First() * binsize_xy.First())
        && ( histogram.GetBoundariesY().Second() == minmax_y.First() + number_of_bins_xy.Second() * binsize_xy.Second())
        && ( histogram.GetNumberOfBinsX() == number_of_bins_xy.First())
        && ( histogram.GetNumberOfBinsY() == number_of_bins_xy.Second()),
        "constructing histogram from bin size xy, number bins xy and starting bin xy did not work"
      );

      //Vector of bin-values (middle point of each bin)
      BCL_MessageStd
      (
        std::string( "these are the binning values - middlepoint of each bin in x and y direction as pair")
      );
      BCL_MessageStd
      (
        util::Format()( histogram.GetBinningXY())
      );

      histogram.CalculateHistogram( values_vector);

      //the core of the histogram without the boundary counts
      BCL_MessageStd
      (
        std::string( "this is just the core  - every count in the bins without the boundaries")
      );
      BCL_MessageStd( util::Format()( histogram.GetHistogram()));

      //test read nd write function
      BCL_MessageStd
      (
        std::string( "this is the test for the read and write function of the Histogram2D class")
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( histogram);
      // read from file
      math::Histogram2D newhistogram;
      ReadBCLObject( newhistogram);
      BCL_MessageStd( util::Format()( newhistogram));

      BCL_MessageStd( std::string( "push back valuepairs 0, 0 and 2, 2 to the histogram"));
      newhistogram.PushBack( storage::VectorND< 2, double>( 0, 0));
      newhistogram.PushBack( storage::VectorND< 2, double>( 2, 2));
      BCL_MessageStd( util::Format()( newhistogram));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathHistogram2D

  const ExampleClass::EnumType ExampleMathHistogram2D::s_Instance
  (
    GetExamples().AddEnum( ExampleMathHistogram2D())
  );

} // namespace bcl

