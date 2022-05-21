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
#include "biol/bcl_biol_blast_profile.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_blast_profile.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolBlastProfile :
    public ExampleInterface
  {
  public:

    ExampleBiolBlastProfile *Clone() const
    {
      return new ExampleBiolBlastProfile( *this);
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
      // initialize ifstream and filenames and vector names
      io::IFStream read;
      io::OFStream write;
      const std::string vector_a_filename( AddExampleInputPathToFilename( e_Biology, "blast_profile_vector_a.txt"));
      const std::string vector_b_filename( AddExampleInputPathToFilename( e_Biology, "blast_profile_vector_b.txt"));
      linal::Vector< double> profile_vector_a( 20), profile_vector_b( 20);
      linal::Vector< double> profile_probabilities_a( 20), profile_probabilities_b( 20);

      BCL_MessageStd( "Reading first profile vector from file");
      // open first vector file
      BCL_ExampleMustOpenInputFile( read, vector_a_filename);
      // read first vector
      read >> profile_vector_a;
      read >> profile_probabilities_a;

      // clear read
      io::File::CloseClearFStream( read);

      BCL_MessageStd( "Reading second profile vector from file");
      // open second vector file
      BCL_ExampleMustOpenInputFile( read, vector_b_filename);
      // read second vector
      read >> profile_vector_b;
      read >> profile_probabilities_b;

      // clear read
      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "Testing default constructor");
      biol::BlastProfile blast_profile_empty;

      // test constructor from a vector of doubles
      BCL_MessageStd( "Testing constructor from a vector");
      biol::BlastProfile blast_profile_a( profile_vector_a);
      biol::BlastProfile blast_profile_b( profile_vector_b);

      // test copy constructor
      BCL_MessageStd( "Testing copy constructor");
      biol::BlastProfile blast_profile_a_copy( blast_profile_a);

      // test clone
      BCL_MessageStd( "Testing clone function");
      util::ShPtr< biol::BlastProfile> sp_blast_profile( blast_profile_a.Clone());

    /////////////////
    // data access //
    /////////////////

      // example check
      BCL_ExampleCheck( blast_profile_empty.GetClassIdentifier(), GetStaticClassName( blast_profile_empty));

      // test the GetProfile( AATYPe) function
      BCL_MessageStd( "Testing GetProfile( AAType) function");
      BCL_ExampleCheck( blast_profile_a.GetProfile( biol::GetAATypes().CYS), profile_vector_a( biol::GetAATypes().CYS));

      // test the GetProfile function
      BCL_MessageStd( "Testing GetProfile function");

      // for blast_profile_a
      BCL_MessageStd( "Profile for blast_profile_a\n" + util::Format()( blast_profile_a.GetProfile()));
      BCL_ExampleCheck( blast_profile_a.GetProfile(), profile_vector_a);

      // for blast_profile_b
      BCL_MessageStd( "Profile for blast_profile_b\n" + util::Format()( blast_profile_b.GetProfile()));
      BCL_ExampleCheck( blast_profile_b.GetProfile(), profile_vector_b);

      // test the GetProfile function
      BCL_MessageStd( "Testing GetProbabilities function");

      // for blast_profile_a
      BCL_MessageStd( "Probabilities for blast_profile_a\n" + util::Format()( blast_profile_a.GetProbabilities()));
      BCL_ExampleIndirectCheckWithinTolerance
      (
        blast_profile_a.GetProbabilities(),
        profile_probabilities_a,
        0.001,
        "Reading BLAST profile"
      );

      // for blast_profile_b
      BCL_MessageStd( "Probabilities for blast_profile_b\n" + util::Format()( blast_profile_b.GetProbabilities()));
      BCL_ExampleIndirectCheckWithinTolerance
      (
        blast_profile_b.GetProbabilities(),
        profile_probabilities_b,
        0.001,
        "Reading BLAST profile"
      );

      // test set_profile function
      BCL_MessageStd( "Testing SetProfile function");
      blast_profile_b.SetProfile( profile_vector_a);
      BCL_ExampleIndirectCheck( blast_profile_b.GetProfile(), blast_profile_a.GetProfile(), "SetProfile");
      BCL_ExampleIndirectCheck( blast_profile_b.GetProbabilities(), blast_profile_a.GetProbabilities(), "SetProfile");

      // compare blast_profile_a and blast_profile_a_copy
      BCL_MessageStd( "compare blast_profile_a and blast_profile_a_copy");
      BCL_ExampleIndirectCheck
      (
        blast_profile_b.GetProfile(),
        blast_profile_a_copy.GetProfile(),
        "Copy constructor"
      );
      BCL_ExampleIndirectCheck
      (
        blast_profile_b.GetProbabilities(),
        blast_profile_a_copy.GetProbabilities(),
        "Copy constructor"
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test CalculateConservation
      const double correct_conservation( 0.264231);
      BCL_ExampleCheckWithinTolerance( blast_profile_a.CalculateConservation(), correct_conservation, 0.000001);

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( blast_profile_a);
      BCL_MessageVrb( "read object");
      biol::BlastProfile blast_profile_read;
      ReadBCLObject( blast_profile_read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_ExampleIndirectCheck( blast_profile_read.GetProfile(), blast_profile_a.GetProfile(), "I/O symmetry");
      BCL_ExampleIndirectCheckWithinTolerance
      (
        blast_profile_a.GetProbabilities(),
        blast_profile_read.GetProbabilities(),
        0.001,
        "I/O symmetry"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleBiolBlastProfile

  const ExampleClass::EnumType ExampleBiolBlastProfile::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolBlastProfile())
  );

} // namespace bcl
