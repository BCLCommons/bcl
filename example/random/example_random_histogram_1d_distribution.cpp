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
#include "random/bcl_random_histogram_1d_distribution.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_random_histogram_1d_distribution.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRandomHistogram1DDistribution :
    public ExampleInterface
  {
  public:

    ExampleRandomHistogram1DDistribution *Clone() const
    {
      return new ExampleRandomHistogram1DDistribution( *this);
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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstantPhiPsi
    //! TODO: add a brief comment
    //! TODO: add an general comment to this class
    //!
    //! @author alexanns
    //! @date Sep 5, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class RandomNumberGenerator :
      public random::DistributionInterface
    {

    private:

    //////////
    // data //
    //////////

      double m_RandomDouble;

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      RandomNumberGenerator( const double RANDOM_DOUBLE) :
        m_RandomDouble( RANDOM_DOUBLE)
      {
      }

      RandomNumberGenerator *Clone() const
      {
        return new RandomNumberGenerator( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name if this function is overwritten
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief sets the seed
      //! @param SEED seed to be used
      //! @return the seed, that is was set too
      uint64_t SetSeed( const uint64_t SEED)
      {
        return 0;
      }

      //! @brief get the seed
      //! @return the seed, that was used to start the rng
      uint64_t GetSeed() const
      {
        return 0;
      }

      //! @brief default range for random double
      //! @return range, in which Double() generates a random number
      const math::Range< double> &GetDoubleRange() const
      {
        static const math::Range< double> s_default_range( math::RangeBorders::e_LeftClosed, 0.0, 1.0, math::RangeBorders::e_RightOpen);
        return s_default_range;
      }

      uint64_t Unsigned64BitInt() const
      {
        return uint64_t( 1);
      }

      //! @brief random double
      //! @return random number in double range
      double Double() const
      {
        return m_RandomDouble;
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class RandomNumberGenerator

    int Run() const
    {
      // create histogram that will be used as the probability distribution
      // has bins of size 1 with centers at 0.5, 1.5, 2.5
      // the lower boundary is 0.0, the upper boundary is 3.0
      math::Histogram histogram( 0, 1, 3);

      // the 0.5 bin has a count of 4
      histogram.PushBack( 0.4);
      histogram.PushBack( 0.6);
      histogram.PushBack( 0.7);
      histogram.PushBack( 0.5);
      // the 1.5 bin has a count of 2
      histogram.PushBack( 1.1);
      histogram.PushBack( 1.3);
      // the 2.5 bin has a count of 1
      histogram.PushBack( 2.4);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      random::Histogram1DDistribution def_constr;

      // constructor taking parameters
      RandomNumberGenerator rng( 0.5);
      random::Histogram1DDistribution param_constr( histogram, rng);
      {
        const size_t histogram_index( param_constr.DetermineRandomCase());
        BCL_MessageDbg( "histogram_index is " + util::Format()( histogram_index));
        BCL_ExampleCheck( histogram_index, 0);
      }

      // test clone constructor
      util::ShPtr< random::Histogram1DDistribution> clone_constr( param_constr.Clone());
      {
        const size_t histogram_index( clone_constr->DetermineRandomCase());
        BCL_MessageDbg( "histogram_index is " + util::Format()( histogram_index));
        BCL_ExampleCheck( histogram_index, 0);
      }

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< random::Histogram1DDistribution>(), clone_constr->GetClassIdentifier());

      // test GetDistributionVector
      const linal::Vector< double> distribution( clone_constr->GetDistributionVector());
      BCL_ExampleCheck( distribution.GetSize(), 3);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test DetermineRandomCase
      {
        RandomNumberGenerator rng( 0.58);
        random::Histogram1DDistribution param_constr( histogram, rng);
        const size_t histogram_index( param_constr.DetermineRandomCase());
        BCL_MessageDbg( "histogram_index is " + util::Format()( histogram_index));
        BCL_ExampleCheck( histogram_index, 1);
      }
      {
        RandomNumberGenerator rng( 0.86);
        random::Histogram1DDistribution param_constr( histogram, rng);
        const size_t histogram_index( param_constr.DetermineRandomCase());
        BCL_MessageDbg( "histogram_index is " + util::Format()( histogram_index));
        BCL_ExampleCheck( histogram_index, 2);
      }

      // test DetermineRandomCase taking boundary and size to directly get bin name
      {
        RandomNumberGenerator rng( 0.80);
        random::Histogram1DDistribution param_constr( histogram, rng);
        const double bin( param_constr.DetermineRandomCase( histogram.GetBoundaries().First(), histogram.GetBinSize()));
        BCL_MessageDbg( "bin is " + util::Format()( bin));
        BCL_ExampleCheck( bin, 1.5);
      }
      // test DetermineRandomCase taking boundary and size to directly get bin name
      {
        RandomNumberGenerator rng( 0.89);
        random::Histogram1DDistribution param_constr( histogram, rng);
        const double bin( param_constr.DetermineRandomCase( histogram.GetBoundaries().First(), histogram.GetBinSize()));
        BCL_MessageDbg( "bin is " + util::Format()( bin));
        BCL_ExampleCheck( bin, 2.5);
      }
      // test DetermineRandomCase taking boundary and size to directly get bin name
      {
        RandomNumberGenerator rng( 0.20);
        random::Histogram1DDistribution param_constr( histogram, rng);
        const double bin( param_constr.DetermineRandomCase( histogram.GetBoundaries().First(), histogram.GetBinSize()));
        BCL_MessageDbg( "bin is " + util::Format()( bin));
        BCL_ExampleCheck( bin, 0.5);
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // test Write and read
      RandomNumberGenerator rng_b( 0.5);
      random::Histogram1DDistribution write_distr( histogram, rng);
      WriteBCLObject( write_distr);
      random::Histogram1DDistribution read_distr;
      ReadBCLObject( read_distr);
      BCL_ExampleCheck( read_distr.GetDistributionVector().GetSize(), 3);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRandomHistogram1DDistribution

  const ExampleClass::EnumType ExampleRandomHistogram1DDistribution::s_Instance
  (
    GetExamples().AddEnum( ExampleRandomHistogram1DDistribution())
  );

} // namespace bcl
