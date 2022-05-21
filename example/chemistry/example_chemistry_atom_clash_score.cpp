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
#include "chemistry/bcl_chemistry_atom_clash_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_atom_clash_score.cpp
  //!
  //! @author kothiwsk
  //! @date Aug 28, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAtomClashScore :
    public ExampleInterface
  {
  public:

    ExampleChemistryAtomClashScore *Clone() const
    {
      return new ExampleChemistryAtomClashScore( *this);
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

      // initialize sdf filename
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "1_3_pentadiene_E.sdf"));
      ensemble.Append( chemistry::FragmentEnsemble( input));
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "azulene.clean.sdf"));
      ensemble.Append( chemistry::FragmentEnsemble( input));
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "taxol.clashed.sdf"));
      ensemble.Append( chemistry::FragmentEnsemble( input));
      io::File::CloseClearFStream( input);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // consider hydrogens, tolerance of 0.1
      chemistry::AtomClashScore score_h_medium_tol( true, true);
      chemistry::AtomClashScore score_h_low_tol( false, true);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      storage::List< chemistry::FragmentComplete>:: const_iterator itr( ensemble.GetMolecules().Begin());
      BCL_ExampleCheckWithinAbsTolerance( score_h_medium_tol( *itr), double( 0.0), 0.01);
      BCL_ExampleCheckWithinAbsTolerance( score_h_low_tol( *itr), double( 0), 0.01);

      // very simple structure
      ++itr;
      BCL_ExampleCheckWithinAbsTolerance( score_h_medium_tol( *itr), double( 0), 0.01);
      BCL_ExampleCheckWithinAbsTolerance( score_h_low_tol( *itr), double( 0), 0.01);

      // azulene is a simple planar aromatic ring. The (mostly hydrogen) clashes are a
      ++itr;
      BCL_ExampleCheckWithinAbsTolerance( score_h_medium_tol( *itr), double( 0.0), 0.01);
      BCL_ExampleCheckWithinAbsTolerance( score_h_low_tol( *itr), double( 0.0), 0.01);

      // check the clashed taxol
      ++itr;
      BCL_ExampleCheckWithinAbsTolerance( score_h_medium_tol( *itr), double( 0.107), 0.01);
      BCL_ExampleCheckWithinAbsTolerance( score_h_low_tol( *itr), double( 0.030), 0.01);

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryAtomClashScore

  const ExampleClass::EnumType ExampleChemistryAtomClashScore::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryAtomClashScore())
  );

} // namespace bcl
