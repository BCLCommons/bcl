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
#include "math/bcl_math_object_probability_distribution.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_object_probability_distribution.cpp
  //!
  //! @author mueller
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathObjectProbabilityDistribution :
    public ExampleInterface
  {
  public:

    ExampleMathObjectProbabilityDistribution *Clone() const
    { return new ExampleMathObjectProbabilityDistribution( *this);}

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
       storage::Pair< double, std::string>::s_Instance.IsDefined();

       // create a SharedPointerVector
       const storage::Vector< typename math::ObjectProbabilityDistribution< descriptor::StringSequence>::Assignment>
       random_distribution
       (
         1,
         math::ObjectProbabilityDistribution< descriptor::StringSequence>::Assignment
         (
           0.7,
           descriptor::StringSequence( "yellow")
         )
       );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor of a ObjectProbabilityDistribution object
      math::ObjectProbabilityDistribution< descriptor::StringSequence> default_distribution;

      // make sure the default object has no latent probability
      BCL_ExampleCheck( math::ObjectProbabilityDistribution< descriptor::StringSequence>().CalculateSum(), 0.0);

      // constructor of a ObjectProbabilityDistribution object from members
      math::ObjectProbabilityDistribution< descriptor::StringSequence> distribution( random_distribution);

      BCL_ExampleIndirectCheck( distribution.CalculateSum(), 0.7, "constructor from members");

      // cloning of a ObjectProbabilityDistribution object
      util::ShPtr< math::ObjectProbabilityDistribution< descriptor::StringSequence> > clone_distribution( distribution.Clone());

      BCL_ExampleIndirectCheck( clone_distribution->CalculateSum(), 0.7, "Clone()");

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      BCL_ExampleCheck
      (
        GetStaticClassName< math::ObjectProbabilityDistribution< descriptor::StringSequence> >(),
        "bcl::math::ObjectProbabilityDistribution<bcl::descriptor::StringSequence>"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< math::ObjectProbabilityDistribution< descriptor::StringSequence> >(),
        clone_distribution->GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

      // push back a new object with its probability as a storage::Pair
      distribution.PushBack( 0.5, descriptor::StringSequence( "blue"));
      distribution.PushBack( 0.3, descriptor::StringSequence( "green"));

      // determine a random case from the distribution and output it
      BCL_ExampleCheck( distribution.DetermineRandomCase().GetString(), "blue");

      // calculate the sum of probabilities and display it
      BCL_ExampleCheck( distribution.CalculateSum(), 1.5);

      // write out object
      WriteBCLObject( distribution);

      // read object
      math::ObjectProbabilityDistribution< descriptor::StringSequence> read_distribution;
      ReadBCLObject( read_distribution);

      BCL_ExampleIndirectCheck( read_distribution.CalculateSum(), 1.5, "I/O");

      // output the distribution object
      BCL_MessageStd( " The entire distribution object is: " + util::Format()( distribution));

      // change the probability of sp_random_object_green and redetermine random case and display it
      distribution.PushBack( 999.7, descriptor::StringSequence( "green"));

      BCL_ExampleIndirectCheck
      (
        distribution.CalculateSum(),
        1001.2,
        "effect of changing probability of green to 1000"
      );

      // determine a random case from the distribution and output it
      BCL_ExampleIndirectCheck
      (
        distribution.DetermineRandomCase().GetString(),
        "green",
        "effect of changing probability of green to 1000"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathObjectProbabilityDistribution

  const ExampleClass::EnumType ExampleMathObjectProbabilityDistribution::s_Instance
  (
    GetExamples().AddEnum( ExampleMathObjectProbabilityDistribution())
  );

} // namespace bcl
