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
#include "nmr/bcl_nmr_spectrum.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

#undef AddAtom

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_nmr_spectrum.cpp
  //!
  //! @author mueller
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleNmrSpectrum :
    public ExampleInterface
  {
  public:

    ExampleNmrSpectrum *Clone() const
    { return new ExampleNmrSpectrum( *this);}

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
      //preparing an input stream from a local file
      const std::string input_file( "nmrshiftdb_first2.sdf");
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, input_file));

      //setting up the ensemble
      chemistry::FragmentEnsemble ensemble( input);
      chemistry::FragmentComplete &first_molecule( *ensemble.Begin());

      //reading data from the input stream
      io::File::CloseClearFStream( input);

      //size of the ensemble
      BCL_MessageStd( "size: " + util::Format()( ensemble.GetSize()));

      // example for output of a spectrum
      const storage::Map< size_t, util::ShPtr< nmr::Spectrum> > spectra
      (
        nmr::Spectrum::GenerateSpectra( first_molecule)
      );
      BCL_ExampleCheck( spectra.GetSize(), 2);

      // first spectrum
      const nmr::Spectrum first_spectrum( *spectra.Begin()->second);
      BCL_ExampleCheck( first_spectrum.GetSpecType(), nmr::Spectrum::e_13C);
      BCL_MessageStd
      (
        "13C spectrum of first molecule: " + util::Format()( first_spectrum)
      );

      // second spectrum
      const nmr::Spectrum second_spectrum( *( ++spectra.Begin())->second);
      BCL_ExampleCheck( second_spectrum.GetSpecType(), nmr::Spectrum::e_1H);
      BCL_MessageStd
      (
        "1H spectrum of first molecule: " + util::Format()( second_spectrum)
      );

      // add the spectrum as second spectrum
      first_spectrum.AddSpectrum( first_molecule);

      // write the molecule
      const std::string nmr_molecule_out( AddExampleOutputPathToFilename( first_spectrum, "two_spectra_out.sdf"));

      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, nmr_molecule_out);
      first_molecule.WriteMDL( output);
      io::File::CloseClearFStream( output);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleNmrSpectrum

  const ExampleClass::EnumType ExampleNmrSpectrum::s_Instance
  (
    GetExamples().AddEnum( ExampleNmrSpectrum())
  );

} // namespace bcl
