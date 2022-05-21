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
#include "math/bcl_math_histogram_3d.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_histogram_3d.cpp
  //!
  //! @author mendenjl
  //! @date Jan 05, 2017
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathHistogram3D :
    public ExampleInterface
  {
  public:

    ExampleMathHistogram3D *Clone() const
    {
      return new ExampleMathHistogram3D( *this);
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
      storage::VectorND< 2, double> minmax_z( 0, 18);
      storage::VectorND< 3, double> binsize_xy( 0.5, 0.25, 2);
      storage::VectorND< 3, size_t> number_of_bins_xy( 6, 12, 18);

      //fill valuesvector with random values

      math::Histogram3D histogram
      (
        storage::VectorND< 3, double>( minmax_x.First(), minmax_y.First(), minmax_z.First()),
        binsize_xy,
        number_of_bins_xy
      );
      for( size_t i( 0); i < number_value_pairs; ++i)
      {
        histogram.PushBack
        (
          random::GetGlobalRandom().Random< double>( minmax_x.First() - 2, minmax_x.Second() + 2),
          random::GetGlobalRandom().Random< double>( minmax_y.First() - 2, minmax_y.Second() + 2),
          random::GetGlobalRandom().Random< double>( minmax_z.First() - 2, minmax_z.Second() + 2)
        );
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //instantiate histogram2D
      BCL_MessageStd( std::string( "building and calculating histogram"));

      BCL_Example_Check
      (
        ( histogram.GetBinSizeXYZ()              == binsize_xy)
        && ( histogram.GetBoundariesX().First()  == minmax_x.First())
        && ( histogram.GetBoundariesY().First()  == minmax_y.First())
        && ( histogram.GetBoundariesZ().First()  == minmax_z.First())
        && ( histogram.GetBoundariesX().Second() == minmax_x.First() + number_of_bins_xy.First() * binsize_xy.First())
        && ( histogram.GetBoundariesY().Second() == minmax_y.First() + number_of_bins_xy.Second() * binsize_xy.Second())
        && ( histogram.GetBoundariesZ().Second() == minmax_z.First() + number_of_bins_xy.Third() * binsize_xy.Third())
        && ( histogram.GetNumberOfBinsX() == number_of_bins_xy.First())
        && ( histogram.GetNumberOfBinsY() == number_of_bins_xy.Second())
        && ( histogram.GetNumberOfBinsZ() == number_of_bins_xy.Third()),
        "constructing histogram from bin size xyz, number bins xyz and starting bin xyz did not work"
      );

      //Vector of bin-values (middle point of each bin)
      BCL_MessageStd
      (
        std::string( "these are the binning values - middlepoint of each bin in x and y direction as pair")
      );
      BCL_MessageStd
      (
        util::Format()( histogram.GetBinningXYZ())
      );

      //the core of the histogram without the boundary counts
      BCL_MessageStd
      (
        std::string( "this is just the core  - every count in the bins without the boundaries")
      );
      BCL_MessageStd( util::Format()( histogram.GetHistogram()));

      math::Histogram3D histogram_b
      (
        storage::VectorND< 3, double>( -0.5, -0.5, -0.5),
        storage::VectorND< 3, double>( 1.0, 1.0, 1.0),
        storage::VectorND< 3, size_t>( 2, 2, 2)
      );
      histogram_b.PushBack( 0, 0, 0, 4);
      histogram_b.PushBack( 0, 0, 1, 5);
      histogram_b.PushBack( 0, 1, 0, 6);
      histogram_b.PushBack( 0, 1, 1, 7);
      histogram_b.PushBack( 1, 0, 0, 8);
      histogram_b.PushBack( 1, 0, 1, 9);
      histogram_b.PushBack( 1, 1, 0, 10);
      histogram_b.PushBack( 1, 1, 1, 11);

      math::Tensor< double> expected_histogram_b
      (
        2,
        2,
        2,
        linal::FillVector< double>( 8, 4, 1).Begin()
      );
      BCL_ExampleIndirectCheck( histogram_b.GetHistogram(), expected_histogram_b, "PushBack");

      BCL_ExampleCheck( histogram_b.Interpolate( 0, 0, 0), 4);
      BCL_ExampleCheck( histogram_b.Interpolate( 0, 0, 1), 5);
      BCL_ExampleCheck( histogram_b.Interpolate( 1, 1, 1), 11);
      BCL_ExampleCheck( histogram_b.Interpolate( 2, 2, 2), 18);
      BCL_ExampleCheck( histogram_b.Interpolate( -1, -1, -1), -3);
      BCL_ExampleCheck( histogram_b.Interpolate( 0.5, 0.5, 0.5), 7.5);
      BCL_ExampleCheck( histogram_b.Interpolate( 0.0, 0.5, 0.5), 5.5);
      BCL_ExampleCheck( histogram_b.Interpolate( 0.25, 0.0, 0.0), 5);

      BCL_ExampleCheck( histogram_b.Value( 0, 0, 0), 4);
      BCL_ExampleCheck( histogram_b.Value( 0, 0, 1), 5);
      BCL_ExampleCheck( histogram_b.Value( 1, 1, 1), 11);
      BCL_ExampleCheck( histogram_b.Value( 2, 2, 2), 11);
      BCL_ExampleCheck( histogram_b.Value( -1, -1, -1), 4);
      BCL_ExampleCheck( histogram_b.Value( 0.5, 0.5, 0.5), 11);
      BCL_ExampleCheck( histogram_b.Value( 0.0, 0.5, 0.5), 7);
      BCL_ExampleCheck( histogram_b.Value( 0.25, 0.0, 0.0), 4);

    //////////////////////
    // input and output //
    //////////////////////

      // Checks Read and Write
      BCL_ExampleIndirectCheck
      (
        ExampleInterface::TestBCLObjectIOForSymmetry( histogram, math::Histogram3D()),
        true,
        "I/O"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathHistogram3D

  const ExampleClass::EnumType ExampleMathHistogram3D::s_Instance
  (
    GetExamples().AddEnum( ExampleMathHistogram3D())
  );

} // namespace bcl

