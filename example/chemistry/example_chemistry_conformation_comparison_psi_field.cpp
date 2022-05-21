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
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "io/bcl_io_file.h"
// external includes

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_conformation_comparison_psi_field.cpp
  //!
  //! @author brownbp1
  //! @date Nov 7, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConformationComparisonPsiField :
    public ExampleInterface
  {
  public:

    ExampleChemistryConformationComparisonPsiField *Clone() const
    { return new ExampleChemistryConformationComparisonPsiField( *this);}

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
      std::string filename_1( AddExampleInputPathToFilename( e_Chemistry, "1dls_moe_version.sdf"));
      std::string filename_2( AddExampleInputPathToFilename( e_Chemistry, "1dhf_moe_version.sdf"));
      io::IFStream input_1( filename_1.c_str());
      io::IFStream input_2( filename_2.c_str());
      BCL_ExampleAssert(input_1.is_open(), true);
      BCL_ExampleAssert(input_2.is_open(), true);

      BCL_MessageStd("Read in 1dls");
      chemistry::FragmentEnsemble mol_1(input_1);
      BCL_MessageStd("Read in 1dhf");
      chemistry::FragmentEnsemble mol_2(input_2);

      //close stream
      io::File::CloseClearFStream( input_1);
      io::File::CloseClearFStream( input_2);

      //set properties and prepare molcules
      descriptor::Combine< chemistry::AtomConformationalInterface, float> properties;
      chemistry::ConformationComparisonPropertyFieldCorrelation property_assigner;
      chemistry::FragmentComplete mol_1_prop(mol_1.GetMolecules().FirstElement());
      chemistry::FragmentComplete mol_2_prop(mol_2.GetMolecules().FirstElement());
      property_assigner.Prepare(mol_1_prop);
      property_assigner.Prepare(mol_2_prop);

      //do alignment and output RMSDX with ConformationComparisonPsiField
      chemistry::ConformationComparisonPsiField comparer;
      double result( comparer( mol_1_prop, mol_2_prop));
      BCL_MessageStd( "PropertyFieldDistance score: " + util::Format()( result));
      BCL_ExampleCheckWithinAbsTolerance
      (
        result,
        0.649482,
        0.05
      );

      //include conformers of mol_1
      std::stringstream error_stream;
      chemistry::SampleConformations sampler;
      sampler.TryRead( util::ObjectDataLabel(), error_stream);
      chemistry::FragmentEnsemble ensemble( sampler( mol_1.GetMolecules().FirstElement()).First());

      // set the keep indices since we are not running operator()()
      storage::Vector< size_t> keep_indices_a, keep_indices_b;

      // a
      for
      (
          auto itr_a( mol_1_prop.GetAtomVector().Begin()), itr_a_end( mol_1_prop.GetAtomVector().End());
          itr_a != itr_a_end;
          ++itr_a
      )
      {
        keep_indices_a.PushBack( mol_1_prop.GetAtomVector().GetAtomIndex( *itr_a));
      }
      comparer.SetKeepIndicesA( keep_indices_a);

      // b
      for
      (
          auto itr_b( mol_2_prop.GetAtomVector().Begin()), itr_b_end( mol_2_prop.GetAtomVector().End());
          itr_b != itr_b_end;
          ++itr_b
      )
      {
        keep_indices_b.PushBack( mol_2_prop.GetAtomVector().GetAtomIndex( *itr_b));
      }
      comparer.SetKeepIndicesB( keep_indices_b);

      // run test
      auto result_ens( comparer.FieldOptimizeOrientation
        (
          mol_1.GetMolecules().FirstElement(),
          mol_2.GetMolecules().FirstElement(),
          400,
          160,
          false,
          0.95,
          ensemble
        ));
      BCL_MessageStd( "PropertyFieldDistance score: " + util::Format()( result_ens.Second()));
      BCL_ExampleCheckWithinAbsTolerance
      (
        result_ens.Second(),
        0.631989,
        0.05
      );

      return 0;
    }

    static const ExampleClass::EnumType ExampleChemistryConformationComparisonPsiField_Instance;

  }; //end ExampleChemistryConformationComparisonPsiField

  const ExampleClass::EnumType ExampleChemistryConformationComparisonPsiField::ExampleChemistryConformationComparisonPsiField_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConformationComparisonPsiField())
  );

} // namespace bcl
