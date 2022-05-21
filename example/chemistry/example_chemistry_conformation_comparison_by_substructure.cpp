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
#include "chemistry/bcl_chemistry_conformation_comparison_by_substructure.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_conformation_comparison_by_substructure.cpp
  //!
  //! @author mendenjl
  //! @date Jun 22, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConformationComparisonBySubstructure :
    public ExampleInterface
  {
  public:

    ExampleChemistryConformationComparisonBySubstructure *Clone() const
    { return new ExampleChemistryConformationComparisonBySubstructure( *this);}

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
      BCL_MessageStd( "All input streams closed");

      chemistry::ConformationComparisonBySubstructure connected_raw_substructure
      (
        false,
        graph::CommonSubgraphIsomorphismBase::e_Connected
      );
      connected_raw_substructure.PrepareEnsemble( ensemble);
      connected_raw_substructure.PrepareEnsemble( napthalene);
      connected_raw_substructure.PrepareEnsemble( taxol);
      connected_raw_substructure.PrepareEnsemble( diazepam);

      // Use first MatchAndRMSD method to find the RMSD between two molecules in the ensemble
      // must move the molecules out of the ensemble and into a ShPtrVector to use MatchAndRMSD on them
      util::SiPtrVector< const chemistry::FragmentComplete> ensemble_vector
      (
        ensemble.GetMolecules().Begin(),
        ensemble.GetMolecules().End()
      );

      const double sub_a( connected_raw_substructure( *ensemble_vector( 1), *ensemble_vector( 0)));
      BCL_MessageStd
      (
        "The Largest common substructure between the first and second elements in the set is: "
        + util::Format()( sub_a)
      );

      BCL_ExampleCheckWithinAbsTolerance
      (
        connected_raw_substructure( *ensemble_vector( 0), *ensemble_vector( 0)),
        ensemble_vector( 0)->GetNumberAtoms(),
        0.5
      );

      // Use first MatchAndRMSD method to find the RMSD between taxol and naphthalene
      BCL_ExampleIndirectCheckWithinTolerance
      (
        connected_raw_substructure
        (
          taxol.GetMolecules().FirstElement(),
          napthalene.GetMolecules().FirstElement()
        ),
        6,
        0.01,
        "ConformationComparisonBySubstructure should be able to compare constitutional-non-isomorphic molecules"
      );

      // calculate the rmsd between molecules 0 and 1 in the vector in both directions to make sure the results are
      // symmetric
      const double sub_1_0( connected_raw_substructure( *ensemble_vector( 1), *ensemble_vector( 0)));
      const double sub_0_1( connected_raw_substructure( *ensemble_vector( 0), *ensemble_vector( 1)));

      // ensure that the MatchAndRMSD function is symmetric
      BCL_ExampleIndirectCheckWithinTolerance( sub_1_0, sub_0_1, 0.001, "Symmetry of MatchAndRMSD");

      return 0;
    }

    static const ExampleClass::EnumType ExampleChemistryConformationComparisonBySubstructure_Instance;

  }; //end ExampleChemistryConformationComparisonBySubstructure

  const ExampleClass::EnumType ExampleChemistryConformationComparisonBySubstructure::ExampleChemistryConformationComparisonBySubstructure_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConformationComparisonBySubstructure())
  );

} // namespace bcl
