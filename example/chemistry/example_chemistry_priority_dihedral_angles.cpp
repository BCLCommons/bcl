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
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_priority_dihedral_angles.cpp
  //!
  //! @author kothiwsk
  //! @date Jul 05, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryPriorityDihedralAngles :
    public ExampleInterface
  {
  public:

    ExampleChemistryPriorityDihedralAngles *Clone() const
    {
      return new ExampleChemistryPriorityDihedralAngles( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

      // load an ensemble directly from an input stream
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));
      chemistry::FragmentEnsemble ensemble( input, sdf::e_Maintain);
      io::File::CloseClearFStream( input);

      storage::Pair< storage::Vector< double>, storage::Vector< storage::VectorND< 4, size_t> > > priority_info
      (
        chemistry::PriorityDihedralAngles()( ensemble.GetMolecules().ExtractElement( 1))
      );

      BCL_ExampleCheck( priority_info.First().GetSize(), size_t( 25));
      BCL_ExampleCheck( math::EqualWithinTolerance( priority_info.First()( 5), double( 89.9986)), true);
      BCL_ExampleCheck( math::EqualWithinTolerance( priority_info.First()( 13), double( 124.999)), true);
      BCL_ExampleCheck( math::EqualWithinTolerance( priority_info.First()( 17), double( 29.7606)), true);

      storage::VectorND< 4, size_t> fifth_vector( 10, 6, 0, 9);
      storage::VectorND< 4, size_t> thirteenth_vector( 20, 15, 3, 2);
      storage::VectorND< 4, size_t> seventeenth_vector( 38, 19, 14, 12);

      BCL_ExampleCheck( priority_info.Second().GetSize(), priority_info.First().GetSize());
      BCL_ExampleCheck( priority_info.Second()( 5), fifth_vector);
      BCL_ExampleCheck( priority_info.Second()( 13), thirteenth_vector);
      BCL_ExampleCheck( priority_info.Second()( 17), seventeenth_vector);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryPriorityDihedralAngles

  const ExampleClass::EnumType ExampleChemistryPriorityDihedralAngles::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryPriorityDihedralAngles())
  );

} // namespace bcl
