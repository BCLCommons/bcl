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
#include "random/bcl_random_histogram_2d_distribution.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_random_histogram_2d_probability_distribution.cpp
  //!
  //! @author rouvelgh
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRandomHistogram2DDistribution :
    public ExampleInterface
  {

  public:

    ExampleRandomHistogram2DDistribution *Clone() const
    {
      return new ExampleRandomHistogram2DDistribution( *this);
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
      // preparation

      // file to be read
      const std::string filename( AddExampleInputPathToFilename( e_Math, "ala_phi_psi.histogram2D"));

      // open input file stream
      io::IFStream read;

      // open IFStream with filename
      BCL_ExampleMustOpenInputFile( read, filename);

      // create histogram2D object
      math::Histogram2D ala_histogram;

      // insert the read file into the histogram object
      read >> ala_histogram;

      // close and clear the  file stream
      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( " Verifying the default constructor")
      random::Histogram2DDistribution default_constructor;

      // initialize ProbabilityDistribution2D object with a histogram2D
      random::Histogram2DDistribution ala_distribution2D( ala_histogram);

      // test constructor from a given file containing a histogram2D for alanine
      BCL_MessageStd( " Verifying the constructor taking a file");

      // test DetermineRandomCase2D
      BCL_MessageStd( "Verifying DetermineRandomCase2D")
      storage::VectorND< 2, size_t> random_case_2d( ala_distribution2D.DetermineRandomCase2D());
      BCL_Example_Check
      (
        util::IsDefined( random_case_2d.First()) &&
        random_case_2d.First() >= 0 &&
        util::IsDefined( random_case_2d.Second()) &&
        random_case_2d.Second() >= 0,
        "The function DetermineRandomCase2D generated unacceptable results of " + util::Format()( random_case_2d.First())
        + " and " + util::Format()( random_case_2d.Second()) + " both of which should always be defined and positive."
      )

      // make a new histogram2D object with an empty histogram2D with size ala_histogram
      math::Histogram2D comparison_histogram( ala_histogram);

      // reset all the counts in the histogram to zero
      comparison_histogram.Reset();

      // fill the empty histogram2D randomly using data from the ala_distribution2D
      const size_t number_iterations( ala_histogram.GetSumOfAllCounts());
      for( size_t counter( 0); counter <= number_iterations; ++counter)
      {
        const storage::VectorND< 2, double> phi_psi
        (
          ala_distribution2D.DetermineRandomCase2D
          (
            comparison_histogram.GetBoundariesX().First(),
            comparison_histogram.GetBoundariesY().First(),
            comparison_histogram.GetBinSizeXY().First(),
            comparison_histogram.GetBinSizeXY().Second()
          )
        );

        comparison_histogram.PushBack( phi_psi);
      }

      // verify the integrity of the generated histogram
      const double comp_rmsd( ( ala_histogram.GetHistogram() - comparison_histogram.GetHistogram()).AsVector().Norm());
      const double max_rmsd( 300.0);
      BCL_Example_Check
      (
        comp_rmsd < max_rmsd,
        "The calculated norm is " + util::Format()( comp_rmsd) + " but should be < "
        + util::Format()( max_rmsd)
      )

    //////////////////////
    // input and output //
    //////////////////////

      // write the comparison histogram to a file
      // file to be written
      const std::string out_name
      (
        AddExampleOutputPathToFilename( comparison_histogram, "comparison_histogram.histogram2D")
      );

      // initialize write stream
      io::OFStream write;

      BCL_ExampleMustOpenOutputFile( write, out_name);
      write << comparison_histogram;
      io::File::CloseClearFStream( write);

      // read object
      random::Histogram2DDistribution read_distribution_2d;
      WriteBCLObject( ala_distribution2D);
      ReadBCLObject( read_distribution_2d);

      // Verify that the read and written objects are identical
      const linal::Matrix< double> difference_matrix( read_distribution_2d.GetMatrix() - ala_distribution2D.GetMatrix());
      const double norm( difference_matrix.AsVector().Norm());
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( 0.0, norm, 1.0E-05),
        "The norm of the difference matrix is: " + util::Format()( norm) + " but should be 0.0"
      )

      // output the read_distribution object
      BCL_MessageDbg( " The entire distribution object is: " + util::Format()( read_distribution_2d));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRandomHistogram2DDistribution

  const ExampleClass::EnumType ExampleRandomHistogram2DDistribution::s_Instance
  (
    GetExamples().AddEnum( ExampleRandomHistogram2DDistribution())
  );

} // namespace bcl

