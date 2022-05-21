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
#include "assemble/bcl_assemble_protein_ensemble.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_protein_ensemble.cpp
  //! TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Feb 11, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleProteinEnsemble :
    public ExampleInterface
  {
    // struct for comparing proteins based on number of residues
    struct SizeCompare
    {
      bool operator()
      (
        const util::ShPtr< assemble::ProteinModel> &MODEL_A, const util::ShPtr< assemble::ProteinModel> &MODEL_B
      )
      {
        return MODEL_A->GetNumberAAs() < MODEL_B->GetNumberAAs();
      }
    }; //< end struct SizeCompare

  public:

    ExampleAssembleProteinEnsemble *Clone() const
    {
      return new ExampleAssembleProteinEnsemble( *this);
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
      // create proteins
      const util::ShPtr< assemble::ProteinModel> ubi
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb")).Clone()
      );
      const util::ShPtr< assemble::ProteinModel> ie
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb")).Clone()
      );
      const std::string lzm_filename( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb"));
      const util::ShPtr< assemble::ProteinModel> lzm( Proteins::GetModel( lzm_filename).Clone());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      assemble::ProteinEnsemble default_constr;
      BCL_ExampleCheck( default_constr.IsEmpty(), true);

      // constructor taking range of iterators
      // create list that can be used
      const util::ShPtrList< assemble::ProteinModel> list( 2, lzm);
      {
        // construct from range of iterators
        const assemble::ProteinEnsemble param_constr( list.Begin(), list.End());
        BCL_ExampleCheck( param_constr.GetSize(), 2);
        BCL_ExampleCheck( *param_constr.Begin(), lzm);
        BCL_ExampleCheck( *++param_constr.Begin(), lzm);
      }

      // clone constructor
      default_constr.InsertElement( ubi);
      util::ShPtr< assemble::ProteinEnsemble> clone_constr( default_constr.Clone());
      BCL_ExampleCheck( clone_constr->GetSize(), 1);

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< assemble::ProteinEnsemble>(), clone_constr->GetClassIdentifier());

      // GetSize
      BCL_ExampleCheck( clone_constr->GetSize(), 1);

      // Begin
      BCL_ExampleCheck( clone_constr->Begin() == --clone_constr->End(), true);

      // End
      BCL_ExampleCheck( ++clone_constr->Begin() == clone_constr->End(), true);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // IsEmpty
      BCL_ExampleCheck( default_constr.IsEmpty(), false);
      BCL_ExampleCheck( assemble::ProteinEnsemble().IsEmpty(), true);

      // InsertElement
      default_constr.InsertElement( ie);
      BCL_ExampleCheck( default_constr.GetSize(), 2);

      // RemoveElement
      default_constr.RemoveElement( default_constr.Begin());
      BCL_ExampleCheck( default_constr.GetSize(), 1);

      // Sort
      clone_constr->InsertElement( ie);
      clone_constr->Sort( SizeCompare());
      BCL_ExampleCheck( clone_constr->GetSize(), 2);
      BCL_ExampleCheck( ( *clone_constr->Begin()), ubi);
      BCL_ExampleCheck( ( *++clone_constr->Begin()), ie);

      // read methods
      storage::Set< sspred::Method> ss_methods( sspred::GetMethods().e_JUFO);
      BCL_ExampleCheck( default_constr.ReadSSPredictions( ss_methods), true);

      // Difference
      { // no difference
        assemble::ProteinEnsemble tmp_ensemble;
        tmp_ensemble.InsertElement( ie);
        tmp_ensemble.InsertElement( ubi);
        const assemble::ProteinEnsemble diff_ensemble( clone_constr->Difference( tmp_ensemble, SizeCompare()));
        BCL_ExampleCheck( diff_ensemble.IsEmpty(), true);
      }
      { // one different in this ensemble
        assemble::ProteinEnsemble tmp_ensemble;
        tmp_ensemble.InsertElement( ie);
        const assemble::ProteinEnsemble diff_ensemble( clone_constr->Difference( tmp_ensemble, SizeCompare()));
        BCL_ExampleCheck( *diff_ensemble.Begin(), ubi);
      }
      { // one different in other ensemble
        assemble::ProteinEnsemble tmp_ensemble;
        tmp_ensemble.InsertElement( ubi);
        tmp_ensemble.InsertElement( ie);
        tmp_ensemble.InsertElement( lzm);
        const assemble::ProteinEnsemble diff_ensemble( clone_constr->Difference( tmp_ensemble, SizeCompare()));
        BCL_ExampleCheck( diff_ensemble.IsEmpty(), true);
      }

      // GetDistanceStatistics
      {
        // make two locators
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_a
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 70)
        );
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_b
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 75)
        );
        const restraint::DataPairwise data_pair( locator_a, locator_b);
        const math::RunningAverageSD< double> mean_sd( clone_constr->GetDistanceStatistics( data_pair));
        BCL_ExampleCheckWithinTolerance( mean_sd.GetAverage(), 13.0283, 0.001);
        BCL_ExampleCheckWithinTolerance( mean_sd.GetStandardDeviation(), 3.20046, 0.001);
      }

      // GetCoordinates
      {
        // make locator
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_a
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 70)
        );
        const storage::Vector< linal::Vector3D> coords( clone_constr->GetCoordinates( *locator_a));
        const storage::Vector< linal::Vector3D> coords_correct
        (
          storage::Vector< linal::Vector3D>::Create
          (
            linal::Vector3D( 29.462, 35.087, 23.43), linal::Vector3D( 5.028, 15.843, 39.253)
          )
        );
        BCL_ExampleCheckWithinTolerance( coords( 0), coords_correct( 0), 0.001);
        BCL_ExampleCheckWithinTolerance( coords( 1), coords_correct( 1), 0.001);
      }

      // GetDistances
      {
        // make two locators
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_a
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 70)
        );
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_b
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 75)
        );
        const restraint::DataPairwise data_pair( locator_a, locator_b);
        const storage::Vector< double> distances( clone_constr->GetDistances( data_pair));
        const storage::Vector< double> distances_correct
        (
          storage::Vector< double>::Create( 16.2287, 9.8278)
        );
        BCL_ExampleCheckWithinTolerance( distances( 0), distances_correct( 0), 0.001);
        BCL_ExampleCheckWithinTolerance( distances( 1), distances_correct( 1), 0.001);
      }

      // GetDistanceChanges
      {
        // make two locators
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_a
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 70)
        );
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_b
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 75)
        );
        const restraint::DataPairwise data_pair( locator_a, locator_b);

        // make other ensemble
        assemble::ProteinEnsemble tmp_ensemble;
        tmp_ensemble.InsertElement( lzm);
        tmp_ensemble.InsertElement( ubi);

        const storage::Vector< double> distance_changes( clone_constr->GetDistanceChanges( data_pair, tmp_ensemble));
        const storage::Vector< double> distance_changes_correct
        (
          storage::Vector< double>::Create( 6.30379, 0, -0.0971222, -6.40091)
        );
        BCL_ExampleCheckWithinTolerance( distance_changes( 0), distance_changes_correct( 0), 0.001);
        BCL_ExampleCheckWithinTolerance( distance_changes( 1), distance_changes_correct( 1), 0.001);
        BCL_ExampleCheckWithinTolerance( distance_changes( 2), distance_changes_correct( 2), 0.001);
        BCL_ExampleCheckWithinTolerance( distance_changes( 3), distance_changes_correct( 3), 0.001);
      }

      // GetDistanceChangesMeanSD one data pair
      {
        // make two locators
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_a
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 70)
        );
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_b
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 75)
        );
        const restraint::DataPairwise data_pair( locator_a, locator_b);

        // make other ensemble
        assemble::ProteinEnsemble tmp_ensemble;
        tmp_ensemble.InsertElement( lzm);
        tmp_ensemble.InsertElement( ubi);

        const math::RunningAverageSD< double> distance_changes( clone_constr->GetDistanceChangesMeanSD( data_pair, tmp_ensemble));
        BCL_ExampleCheckWithinTolerance( distance_changes.GetAverage(), -0.0485611, 0.001);
        BCL_ExampleCheckWithinTolerance( distance_changes.GetStandardDeviation(), 4.49192, 0.001);
      }

      // GetDistanceChangesMeanSD data set
      {
        // make two locators
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_a
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 70)
        );
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_b
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 75)
        );
        const restraint::DataPairwise data_pair_a( locator_a, locator_b);
        // make two more locators
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_c
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 60)
        );
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_d
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 65)
        );
        const restraint::DataPairwise data_pair_b( locator_c, locator_d);

        restraint::DataSetPairwise data_set;
        data_set.Insert( data_pair_a);
        data_set.Insert( data_pair_b);

        // make other ensemble
        assemble::ProteinEnsemble tmp_ensemble;
        tmp_ensemble.InsertElement( lzm);
        tmp_ensemble.InsertElement( ubi);

        const std::multimap< math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean>
          distance_changes( clone_constr->GetDistanceChangesMeanSD( data_set, tmp_ensemble));
        BCL_ExampleCheckWithinTolerance( distance_changes.begin()->first.GetAverage(), -0.0485611, 0.001);
        BCL_ExampleCheckWithinTolerance( distance_changes.begin()->first.GetStandardDeviation(), 4.49192, 0.001);
        BCL_ExampleCheckWithinTolerance( ( --distance_changes.end())->first.GetAverage(), 0.113967, 0.001);
        BCL_ExampleCheckWithinTolerance( ( --distance_changes.end())->first.GetStandardDeviation(), 0.396115, 0.001);
        BCL_ExampleCheck( distance_changes.size(), 2);
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( default_constr);

      // read the object back in
      assemble::ProteinEnsemble read;
      ReadBCLObject( read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      {
        BCL_ExampleCheck( read.GetSize(), 1);
      }

    //////////////////////
    // helper functions //
    //////////////////////

      // GetEnsembleFromFile
      {
        const std::string ensemble_file( AddExampleInputPathToFilename( e_Biology, "2LZM_ensemble.ls"));
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, ensemble_file);
        write << lzm_filename;
        io::File::CloseClearFStream( write);
        const assemble::ProteinEnsemble ensemble( assemble::ProteinEnsemble::GetEnsembleFromFile( ensemble_file, 0, biol::GetAAClasses().e_AAComplete));

        BCL_ExampleCheck( ensemble.GetSize(), 1);
        BCL_ExampleCheck( ensemble.Begin()->IsDefined(), true);
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleProteinEnsemble

  const ExampleClass::EnumType ExampleAssembleProteinEnsemble::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleProteinEnsemble())
  );

} // namespace bcl

