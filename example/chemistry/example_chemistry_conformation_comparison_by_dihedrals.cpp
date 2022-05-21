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
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedrals.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_conformation_comparison_by_dihedrals.cpp
  //!
  //! @author kothiwsk, mendenjl
  //! @date Mar 05, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConformationComparisonByDihedrals :
    public ExampleInterface
  {
  public:

    ExampleChemistryConformationComparisonByDihedrals *Clone() const
    { return new ExampleChemistryConformationComparisonByDihedrals( *this);}

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
      const std::string filename_d( AddExampleInputPathToFilename( e_Chemistry, "corina_diazepam.sdf"));
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
      chemistry::FragmentEnsemble napthalene( input_b, sdf::e_Saturate);
      BCL_MessageStd( "Read in taxol");
      chemistry::FragmentEnsemble taxol( input_c);
      BCL_MessageStd( "Read in diazepam");
      chemistry::FragmentEnsemble diazepam( input_d); // this seems to be the problem

      // close and clear input streams
      io::File::CloseClearFStream( input_a);
      io::File::CloseClearFStream( input_b);
      io::File::CloseClearFStream( input_c);
      io::File::CloseClearFStream( input_d);
      BCL_MessageStd( "All input streams closed");

      // Use first MatchAndRMSD method to find the RMSD between two molecules in the ensemble
      // must move the molecules out of the ensemble and into a ShPtrVector to use MatchAndRMSD on them
      storage::Vector< chemistry::FragmentComplete> ensemble_vector
      (
        ensemble.GetMolecules().Begin(),
        ensemble.GetMolecules().End()
      );

      const double max_difference_dihedral_angles_different_conformations
      (
        chemistry::ConformationComparisonByDihedrals()( ensemble_vector( 1), ensemble_vector( 0))
      );
      BCL_MessageStd
      (
        "The max difference in dihedral angles between the first elements in the set and the second is: "
        + util::Format()( max_difference_dihedral_angles_different_conformations)
      );

      BCL_MessageStd( "dihedral_max_difference_self " + util::Format()( ensemble_vector( 0).GetBondInfo().FirstElement().GetConfigurationalBondType() == ensemble_vector( 0).GetBondInfo().FirstElement().GetConfigurationalBondType()));
      const double dihedral_max_difference_self
      (
        chemistry::ConformationComparisonByDihedrals()( ensemble_vector( 0), ensemble_vector( 0))
      );
      BCL_MessageStd
      (
        "The max difference in dihedral angles between the first elements in the set and itself is: "
        + util::Format()( dihedral_max_difference_self)
      );

      BCL_ExampleIndirectCheckWithinTolerance
      (
        dihedral_max_difference_self, 0.0, 0.001,
        " difference between identical conformations should always be zero"
      );

      // Use first MatchAndRMSD method to find the RMSD between taxol and naphthalene
      BCL_ExampleIndirectCheck
      (
        util::IsDefined
        (
          chemistry::ConformationComparisonByDihedrals()
          (
            taxol.GetMolecules().FirstElement(),
            napthalene.GetMolecules().FirstElement()
          )
        ),
        false,
        "ConformationComparisonByDihedrals should not be able to compare constitutional-non-isomorphic molecules"
      );
      // calculate the rmsd between molecules 0 and 1 in the vector in both directions to make sure the results are
      // symmetric
      const double rmsd_1_0( chemistry::ConformationComparisonByDihedrals()( ensemble_vector( 1), ensemble_vector( 0)));
      const double rmsd_0_1( chemistry::ConformationComparisonByDihedrals()( ensemble_vector( 0), ensemble_vector( 1)));

      // ensure that the MatchAndRMSD function is symmetric
      BCL_ExampleIndirectCheckWithinTolerance( rmsd_1_0, rmsd_0_1, 0.001, "Symmetry of MatchAndRMSD");

      return 0;
    }

    static const ExampleClass::EnumType ExampleChemistryConformationComparisonByDihedrals_Instance;

  }; //end ExampleChemistryConformationComparisonByDihedrals

  const ExampleClass::EnumType ExampleChemistryConformationComparisonByDihedrals::ExampleChemistryConformationComparisonByDihedrals_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConformationComparisonByDihedrals())
  );

} // namespace bcl
