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
#include "random/bcl_random_uniform_distribution.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_random_uniform_distribution.cpp
  //!
  //! @author woetzen, mendenjl
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRandomUniformDistribution :
    public ExampleInterface
  {
  public:

    ExampleRandomUniformDistribution *Clone() const
    { return new ExampleRandomUniformDistribution( *this);}

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
      const size_t seed( 123);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create default random number generator
      random::UniformDistribution rng_default;

      // create random number generator with defined seed
      random::UniformDistribution rng( seed);

      // clone
      util::ShPtr< random::DistributionInterface> ptr( rng.Clone());

    /////////////////
    // data access //
    /////////////////

      storage::Vector< double> poisson_table( 1000);

      for( size_t i( 0); i < 1000; ++i)
      {
        poisson_table( i) = rng.RandomPoisson( 8);
      }

      BCL_Message( util::Message::e_Standard, "Poisson Values: " + util::Format()( poisson_table));

      BCL_ExampleCheck( GetStaticClassName< random::UniformDistribution>(), ptr->GetClassIdentifier());

      // get seed
      BCL_MessageStd( "seed: " + util::Format()( rng.GetSeed()));
      BCL_ExampleIndirectCheck( rng.GetSeed(), seed, "constructor from seed");

      // get range for double random numbers
      BCL_MessageStd( "this is the range in which double random numbers are generated: " + rng.GetDoubleRange().GetString());
      BCL_ExampleCheck( rng.GetDoubleRange(), random::DistributionInterface::GetDefaultDoubleRange());

      // get default seed
      BCL_MessageStd( "default seed: " + util::Format()( random::UniformDistribution::GetDefaultSeed()));

      // define the seed
      BCL_Example_Check
      (
        rng_default.GetSeed() != seed,
        "highly improbable that default constructor was setting seed to: " + util::Format()( seed)
      );
      rng_default.SetSeed( seed);
      BCL_ExampleIndirectCheck( rng_default.GetSeed(), seed, "rng_default.SetSeed( seed)");

      // define seed from commandline flag, which should be equal to the default seed if none was passed over command
      // line
      rng_default.SetSeedFromCommandlineFlag();
      BCL_Example_Check
      (
        rng_default.GetSeed() == random::UniformDistribution::GetDefaultSeed(),
        "the seed was changed in commandline or an error happened: " + util::Format()( rng_default.GetSeed()) + " != " +
        util::Format()( random::UniformDistribution::GetDefaultSeed())
      );

    ///////////////
    // operators //
    ///////////////

      rng_default.SetSeed( seed);

      // assignment
      random::UniformDistribution rng_assigned;
      rng_assigned = rng_default;
      BCL_Example_Check
      (
        rng_default.GetSeed() == rng_assigned.GetSeed() &&
        rng_default.Unsigned64BitInt() == rng_assigned.Unsigned64BitInt(),
        "assigned rng behaves differently the source rng! " + util::Format()( rng_assigned) + "\n" +
        util::Format()( rng_default)
      );
      rng_default.SetSeedFromCommandlineFlag();

    ////////////////
    // operations //
    ////////////////

      // the actual random number, all other random functions, are using this and scale it to their needs
      // note: all messages that print out a random number have to be critical, otherwise other example checks
      // may fail at message level critical
      BCL_MessageCrt( "Random 64 bit unsigned int : " + util::Format()( rng.Unsigned64BitInt()));

      // random numbers in default ranges
      BCL_MessageCrt( "Random double " + rng.GetDoubleRange().GetString()  + ": " + util::Format()( rng.Double()));
      BCL_MessageCrt( "Random float  " + rng.GetFloatRange().GetString()   + ": " + util::Format()( rng.Float()));
      BCL_MessageCrt( "Random 64 bit unsigned int " + rng.GetRange< uint64_t>().GetString() + ": " + util::Format()( rng.Unsigned64BitInt()));
      BCL_MessageCrt( "Random bool   " + rng.GetBooleanRange().GetString() + ": " + util::Format()( rng.Boolean()));
      BCL_MessageCrt( "Random sign : " + util::Format()( rng.Sign()));

      // create different ranges
      const math::Range< double> dcc( math::RangeBorders::e_LeftClosed, 0.0, 1.0, math::RangeBorders::e_RightClosed);
      const math::Range< double> dco( math::RangeBorders::e_LeftClosed, 0.0, 1.0, math::RangeBorders::e_RightOpen);
      const math::Range< double> doo( math::RangeBorders::e_LeftOpen,   0.0, 1.0, math::RangeBorders::e_RightOpen);
      const math::Range<  float> fcc( math::RangeBorders::e_LeftClosed, 0.0, 1.0, math::RangeBorders::e_RightClosed);
      const math::Range<  float> fco( math::RangeBorders::e_LeftClosed, 0.0, 1.0, math::RangeBorders::e_RightOpen);
      const math::Range<  float> foo( math::RangeBorders::e_LeftOpen,   0.0, 1.0, math::RangeBorders::e_RightOpen);
      const math::Range< size_t> s1_range( math::RangeBorders::e_LeftClosed,    0  , 250  , math::RangeBorders::e_RightOpen);
      const math::Range< size_t> s2_range( math::RangeBorders::e_LeftClosed,  250  , 500  , math::RangeBorders::e_RightOpen);
      const math::Range< int>    i1_range( math::RangeBorders::e_LeftClosed, -125  , 125  , math::RangeBorders::e_RightOpen);
      const math::Range< double> d1_range( math::RangeBorders::e_LeftClosed,    0.0, 250.0, math::RangeBorders::e_RightOpen);
      const math::Range< double> d2_range( math::RangeBorders::e_LeftClosed, -125.0, 125.0, math::RangeBorders::e_RightOpen);

      // random numbers in ranges
      BCL_MessageCrt( "random numbers in different ranges");
      BCL_MessageCrt( "Random double " + dcc.GetString() + ": " + util::Format()( rng.Double( dcc)));
      BCL_MessageCrt( "Random double " + dco.GetString() + ": " + util::Format()( rng.Double( dco)));
      BCL_MessageCrt( "Random double " + doo.GetString() + ": " + util::Format()( rng.Double( doo)));
      BCL_MessageCrt( "Random  float " + fcc.GetString() + ": " + util::Format()( rng.Float( fcc)));
      BCL_MessageCrt( "Random  float " + fco.GetString() + ": " + util::Format()( rng.Float( fco)));
      BCL_MessageCrt( "Random  float " + foo.GetString() + ": " + util::Format()( rng.Float( foo)));
      BCL_MessageCrt( "Random size_t " + s1_range.GetString() + ": " + util::Format()( rng.Random( s1_range)));
      BCL_MessageCrt( "Random size_t " + s2_range.GetString() + ": " + util::Format()( rng.Random( s2_range)));
      BCL_MessageCrt( "Random int    " + i1_range.GetString() + ": " + util::Format()( rng.Integer( i1_range)));
      BCL_MessageCrt( "Random double " + d1_range.GetString() + ": " + util::Format()( rng.Double( d1_range)));
      BCL_MessageCrt( "Random double " + d2_range.GetString() + ": " + util::Format()( rng.Double( d2_range)));

      // create sequence of random numbers
      const math::Range< size_t> sequence_range( 0, 9);
      std::string random_sequence;
      for( size_t i( 0); i < 25; ++i)
      {
        random_sequence += util::Format()( rng.Random( sequence_range)) + " ";
      }
      BCL_MessageStd( "Random size_ts " + sequence_range.GetString() + ": " + random_sequence);

      // histograms in order to check the distribution of the random values
      math::Histogram hist_s1(    0, 25, 10);
      math::Histogram hist_s2(  250, 25, 10);
      math::Histogram hist_i1( -125, 25, 10);
      math::Histogram hist_d1(    0, 25, 10);
      math::Histogram hist_d2( -125, 25, 10);

      //check that the ranges are not violated and the random numbers are well distributed
      bool s1_within( true);
      bool s2_within( true);
      bool i1_within( true);
      bool d1_within( true);
      bool d2_within( true);

      // fill histograms by generate lots of random numbers
      for( size_t i( 0); i < 100000; ++i)
      {
        // generate random numbers
        const size_t s1( rng.Random( s1_range));
        hist_s1.PushBack( s1);
        const size_t s2( rng.Random( s2_range));
        hist_s2.PushBack( s2);
        const int    i1( rng.Integer( i1_range));
        hist_i1.PushBack( i1);
        const double d1( rng.Double( d1_range));
        hist_d1.PushBack( d1);
        const double d2( rng.Double( d2_range));
        hist_d2.PushBack( d2);

        // check that umbers are within given interval
        s1_within &= s1_range.IsWithin( s1);
        s2_within &= s2_range.IsWithin( s2);
        i1_within &= i1_range.IsWithin( i1);
        d1_within &= d1_range.IsWithin( d1);
        d2_within &= d2_range.IsWithin( d2);
      }

      BCL_Example_Check( s1_within, "random size_t not in range " + s1_range.GetString());
      BCL_Example_Check( s2_within, "random size_t not in range " + s2_range.GetString());
      BCL_Example_Check( i1_within, "random int not in range " + i1_range.GetString());
      BCL_Example_Check( d1_within, "random double not in range " + d1_range.GetString());
      BCL_Example_Check( d2_within, "random double not in range " + d2_range.GetString());

      BCL_MessageStd( "random size_t distribution in range " + s1_range.GetString());
      BCL_MessageStd( util::Format()( hist_s1));

      BCL_MessageStd( "random size_t distribution in range " + s2_range.GetString());
      BCL_MessageStd( util::Format()( hist_s2));

      BCL_MessageStd( "random    int distribution in range " + i1_range.GetString());
      BCL_MessageStd( util::Format()( hist_i1));

      BCL_MessageStd( "random double distribution in range " + d1_range.GetString());
      BCL_MessageStd( util::Format()( hist_d1));

      BCL_MessageStd( "random double distribution in range " + d2_range.GetString());
      BCL_MessageStd( util::Format()( hist_d2));

    //////////////////////
    // input and output //
    //////////////////////

      // write the random number generator out
      std::stringstream random_number_generator_out;
      rng.SetSeed( seed);
      io::Serialize::Write( rng, random_number_generator_out);
      BCL_MessageVrb( "rng with seed set to " + util::Format()( seed) + "\n" + util::Format()( rng));

      // ensure that the output is unchanged
      const std::string correct_rng_filename( AddExampleInputPathToFilename( e_Random, "bcl_random_UniformDistribution.bcl"));
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, correct_rng_filename);

      // test that the streams matched
      BCL_ExampleCheck( io::File::StreamsMatch( random_number_generator_out, input), true);

      // close the correct file
      io::File::CloseClearFStream( input);

      //build up buckets to check if the random number generator produces well distributed numbers
      BCL_MessageStd
      (
        "build up buckets to check if the random number generator produces well distributed numbers"
      );
      std::map< size_t, size_t> bucket_orig;
      for( size_t i( 0); i < 1000000; ++i)
      {
        bucket_orig[ rng.Random( math::Range< size_t>( 0, 10))]++;
      }

      std::string buckets_string;
      for( size_t i( 0); i < 10; ++i)
      {
        buckets_string += util::Format()( bucket_orig[ i]) + " ";
      }
      BCL_MessageStd( buckets_string);

      //make a new random seed, to see that the random numbers will change
      BCL_MessageStd( "set a new random seed, to see that the random numbers will change");
      rng.SetSeed( 1234);
      bucket_orig.clear();
      for( size_t i( 0); i < 1000000; ++i)
      {
        bucket_orig[ rng.Random( math::Range< size_t>( 0, 10))]++;
      }

      buckets_string = "";
      for( size_t i( 0); i < 10; ++i)
      {
        buckets_string += util::Format()( bucket_orig[ i]) + " ";
      }
      BCL_MessageStd( buckets_string);

      // Finally, write out many random numbers of different types
      // The bcl guarantees that with the same seed, we'll always get the same set of random numbers
      std::stringstream random_numbers;
      rng.SetSeed( 1234);
      random_numbers << '\n' << "64 bit unsigned ints: ";
      for( size_t i( 0); i < 128; ++i)
      {
        random_numbers << rng.Unsigned64BitInt() << ' ';
      }
      random_numbers << '\n' << "Doubles: ";
      for( size_t i( 0); i < 128; ++i)
      {
        random_numbers << rng.Double() << ' ';
      }
      random_numbers << '\n' << "Floats: ";
      for( size_t i( 0); i < 128; ++i)
      {
        random_numbers << rng.Float() << ' ';
      }
      random_numbers << '\n' << "Bools: ";
      for( size_t i( 0); i < 128; ++i)
      {
        random_numbers << rng.Boolean() << ' ';
      }
      random_numbers << '\n' << "Signs: ";
      for( size_t i( 0); i < 128; ++i)
      {
        random_numbers << rng.Sign() << ' ';
      }
      random_numbers << '\n';

      BCL_MessageStd
      (
        "Random numbers of several types:\n"
        + random_numbers.str()
      );

      // check the output against the random numbers in the example_random_numbers.txt filename
      const std::string example_random_numbers_filename
      (
        AddExampleInputPathToFilename( e_Random, "example_random_numbers.txt")
      );
      BCL_ExampleMustOpenInputFile( input, example_random_numbers_filename);
      BCL_ExampleIndirectCheck
      (
        io::File::StreamsMatch( input, random_numbers), true,
        "random number consistency on this platform"
      );

      BCL_MessageStd( "Random numbers output: " + random_numbers.str());
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRandomUniformDistribution

  const ExampleClass::EnumType ExampleRandomUniformDistribution::s_Instance
  (
    GetExamples().AddEnum( ExampleRandomUniformDistribution())
  );

} // namespace bcl
