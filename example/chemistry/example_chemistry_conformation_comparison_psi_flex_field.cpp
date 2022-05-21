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
#include "chemistry/bcl_chemistry_conformation_comparison_psi_flex_field.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_set.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_implementation.h"

// external includes

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_conformation_comparison_psi_flex_field.cpp
  //!
  //! @author brownbp1
  //! @date Nov 7, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConformationComparisonPsiFlexField :
    public ExampleInterface
  {
  public:

    ExampleChemistryConformationComparisonPsiFlexField *Clone() const
    { return new ExampleChemistryConformationComparisonPsiFlexField( *this);}

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
      std::string filename_1( AddExampleInputPathToFilename( e_Chemistry, "1dls_moe_version_corina.sdf"));
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

      std::stringstream error_stream;
      chemistry::SampleConformations sampler;
      sampler.TryRead( util::ObjectDataLabel(), error_stream);

      //set properties and prepare molcules
      descriptor::Combine< chemistry::AtomConformationalInterface, float> properties;
      chemistry::ConformationComparisonPropertyFieldCorrelation property_assigner;
      chemistry::FragmentComplete mol_1_prop(mol_1.GetMolecules().FirstElement());
      chemistry::FragmentComplete mol_2_prop(mol_2.GetMolecules().FirstElement());
      property_assigner.Prepare(mol_1_prop);
      property_assigner.Prepare(mol_2_prop);

      chemistry::FragmentEnsemble ensemble_1( sampler(mol_1_prop).First());
      chemistry::FragmentEnsemble ensemble_2( sampler(mol_2_prop).First());

      //do alignment and output RMSDX with ConformationComparisonPsiField
      chemistry::ConformationComparisonPsiFlexField comparer;

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

      storage::Vector< storage::Triplet< chemistry::FragmentComplete, chemistry::FragmentComplete, double> >
      result_semiflex
      (
        comparer.FieldOptimizeOrientationFlex
        (
          ensemble_1,
          ensemble_2,
          50,
          200,
          50,
          50,
          20,
          50,
          20,
          0.75,
          0.50
        )
      );
      BCL_MessageStd( "Semi-flexible alignment PropertyFieldDistance score: " + util::Format()( result_semiflex( 0).Third()));
      BCL_ExampleCheck
      (
        result_semiflex( 0).Third() <= 0.8,
        true
      );
      return 0;
    }

    static const ExampleClass::EnumType ExampleChemistryConformationComparisonPsiFlexField_Instance;

  }; //end ExampleChemistryConformationComparisonPsiFlexField

  const ExampleClass::EnumType ExampleChemistryConformationComparisonPsiFlexField::ExampleChemistryConformationComparisonPsiFlexField_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConformationComparisonPsiFlexField())
  );

} // namespace bcl
