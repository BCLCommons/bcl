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
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_conformation_comparison_by_dihedral_bins.cpp
  //!
  //! @author kothiwsk
  //! @date Nov 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConformationComparisonByDihedralBins :
    public ExampleInterface
  {
  public:

    ExampleChemistryConformationComparisonByDihedralBins *Clone() const
    { return new ExampleChemistryConformationComparisonByDihedralBins( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      //preparing input streams
      const std::string filename_a( AddExampleInputPathToFilename( e_Chemistry, "CSD_first10sixMemberedRings.sdf"));
      const std::string filename_b( AddExampleInputPathToFilename( e_Chemistry, "naphthalene.sdf"));
      const std::string filename_c( AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      const std::string filename_d( AddExampleInputPathToFilename( e_Chemistry, "pentadiene_conformers.sdf"));
      io::IFStream input_a( filename_a.c_str());
      io::IFStream input_b( filename_b.c_str());
      io::IFStream input_c( filename_c.c_str());
      io::IFStream input_d( filename_d.c_str());

      BCL_ExampleAssert( input_a.is_open(), true);
      BCL_ExampleAssert( input_b.is_open(), true);
      BCL_ExampleAssert( input_c.is_open(), true);
      BCL_ExampleAssert( input_d.is_open(), true);

      //read data from the input streams and store
      BCL_MessageStd( "Read in the ensemble of 10 molecules");
      chemistry::FragmentEnsemble ensemble( input_a);
      BCL_MessageStd( "Read in naphthalene");
      chemistry::FragmentEnsemble napthalene( input_b);
      BCL_MessageStd( "Read in taxol");
      chemistry::FragmentEnsemble taxol( input_c);
      BCL_MessageStd( "Read in diazepam");
      chemistry::FragmentEnsemble diazepam( input_d); // this seems to be the problem

      // close and clear input streams
      io::File::CloseClearFStream( input_a);
      io::File::CloseClearFStream( input_b);
      io::File::CloseClearFStream( input_c);
      io::File::CloseClearFStream( input_d);

      // compare two same molecules
      const double dihedral_max_difference_self
      (
        chemistry::ConformationComparisonByDihedralBins( 30.0)
        (
          ensemble.GetMolecules().FirstElement(), ensemble.GetMolecules().FirstElement()
        )
      );

      BCL_ExampleIndirectCheck
      (
        dihedral_max_difference_self,
        0,
        " difference between identical conformations should always be zero"
      );

      // compare two different molecules
      BCL_ExampleIndirectCheck
      (
        util::IsDefined
        (
          chemistry::ConformationComparisonByDihedralBins( 30.0)
          (
            taxol.GetMolecules().FirstElement(),
            napthalene.GetMolecules().FirstElement()
          )
        ),
        0,
        "CompareConformationsByDihedralBins should not be able to compare constitutional-non-isomorphic molecules"
      );

      // check for symmetry
      const double compare_1_0( chemistry::ConformationComparisonByDihedralBins( 30.0)( ensemble.GetMolecules().FirstElement(), ensemble.GetMolecules().LastElement()));
      const double compare_0_1( chemistry::ConformationComparisonByDihedralBins( 30.0)( ensemble.GetMolecules().LastElement(), ensemble.GetMolecules().FirstElement()));

      BCL_MessageStd( " compare_1_0 : " + util::Format()( compare_0_1));
      BCL_MessageStd( " compare_1_0 : " + util::Format()( compare_1_0));
      // ensure that the  function is symmetric
      BCL_ExampleIndirectCheck
      (
        compare_1_0,
        compare_0_1,
        "Symmetry of comparison"
      );

//    Test the DetermineDihedralKey function
      chemistry::ConformationComparisonByDihedralBins comparer( double( 30));
      chemistry::FragmentComplete diazepam_mol( diazepam.GetMolecules().FirstElement());
      diazepam_mol.RemoveH();
      linal::Vector< double> abs_dihedrals( diazepam_mol.GetDihedralAngles());
      storage::Vector< size_t> keys_diazepam( comparer.DetermineDihedralKeys( abs_dihedrals, false, false));

      storage::Vector< size_t> dihedral_keys_correct( storage::Vector< size_t>::Create( 6, 6));

      BCL_ExampleCheck( keys_diazepam, dihedral_keys_correct);

      storage::Vector< size_t> bins( storage::Vector< size_t>::Create( 1, 2, 3, 4, 5, 6, 7, 8, 11, 12));

      storage::Vector< storage::Pair< double, double> > angle_bounds;

      for
      (
        storage::Vector< size_t>::const_iterator itr( bins.Begin()), itr_end( bins.End());
        itr != itr_end;
        ++itr
      )
      {
        angle_bounds.PushBack( comparer.KeyToAngle( *itr));
      }

      storage::Vector< storage::Pair< double, double> > angle_bounds_correct
      (
        storage::Vector< storage::Pair< double, double> >::Create
        (
          storage::Pair< double, double>( 15, 45),
          storage::Pair< double, double>( 45, 75),
          storage::Pair< double, double>( 75, 105),
          storage::Pair< double, double>( 105, 135),
          storage::Pair< double, double>( 135, 165),
          storage::Pair< double, double>( 165, -165),
          storage::Pair< double, double>( -165, -135),
          storage::Pair< double, double>( -135, -105),
          storage::Pair< double, double>( -45, -15),
          storage::Pair< double, double>( -15, 15)
        )
      );

      BCL_ExampleCheck( angle_bounds, angle_bounds_correct);
      return 0;
    }

    static const ExampleClass::EnumType ExampleChemistryCompareConformationsByDihedralBins_Instance;

  }; //end ExampleChemistryConformationComparisonByDihedralBins

  const ExampleClass::EnumType ExampleChemistryConformationComparisonByDihedralBins::ExampleChemistryCompareConformationsByDihedralBins_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConformationComparisonByDihedralBins())
  );

} // namespace bcl
