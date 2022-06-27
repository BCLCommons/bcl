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
#include "chemistry/bcl_chemistry_mol_align_by_parts.h"

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
  //! @example example_chemistry_mol_align_by_parts.cpp
  //!
  //! @author ben
  //! @date Aug 15, 2021
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryMolAlignByParts :
    public ExampleInterface
  {
  public:

    ExampleChemistryMolAlignByParts *Clone() const
    { return new ExampleChemistryMolAlignByParts( *this);}

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
//      std::string filename_1( "/hd1/brownbp1/MCR_test/her2/h/reaction_products.sdf");
//      std::string filename_1( "/home/ben/Projects/bcl_tests/molalign_by_parts/her2/h/reaction_products.sdf");
//      std::string filename_1( "/hd1/brownbp1/MCR_test/her2/h/reaction_products.screen.unique.sdf");
//      std::string filename_1( "/hd1/brownbp1/MCR_test/her2/h/lig.corina.sdf");
//      std::string filename_2( "/hd1/brownbp1/MCR_test/her2/h/lig.sdf");
//      std::string filename_2( "/home/ben/Projects/bcl_tests/molalign_by_parts/her2/h/lig.sdf");
//      std::string filename_2( "/home/ben/Projects/bcl_tests/molalign_by_parts/her2/h/lig.rigid_split.sdf");
//      std::string filename_2( "/home/ben/Projects/bcl_tests/molalign_by_parts/her2/h/ligs.seq_0.sdf");
//      std::string filename_2( "/home/ben/Projects/bcl_tests/molalign_by_parts/her2/h/ligs.seq_1.sdf");
//      std::string filename_2( "/hd1/brownbp1/MCR_test/her2/h/lig.rigid_split.sdf");
//      io::IFStream input_1( filename_1.c_str());
//      io::IFStream input_2( filename_2.c_str());
//      BCL_ExampleAssert(input_1.is_open(), true);
//      BCL_ExampleAssert(input_2.is_open(), true);
//
//      BCL_MessageStd("Read in " + util::Format()( filename_1));
//      chemistry::FragmentEnsemble mol_1(input_1);
//      BCL_MessageStd("Read in " + util::Format()( filename_2));
//      chemistry::FragmentEnsemble mol_2(input_2);
//
//      //close stream
//      io::File::CloseClearFStream( input_1);
//      io::File::CloseClearFStream( input_2);
//
//      // setup our align by parts object
//      chemistry::MolAlignByParts aligner;
//      storage::Vector< size_t> fixed_indices
//      (
//        util::SplitStringToNumerical< size_t>( util::TrimString( "15,14,13,16,11,12,2,49,50,48,45,46,47,40,41"), " \t\n\r,")
//      );
//      aligner.SetFixedAtoms(fixed_indices);
//
//      //do alignment and output score
//      for
//      (
//          auto mol_itr( mol_1.GetMolecules().Begin()), mol_itr_end( mol_1.GetMolecules().End());
//          mol_itr != mol_itr_end;
//          ++mol_itr
//      )
//      {
//        auto result( aligner.Align(*mol_itr, mol_2));
//      }
//      BCL_ExampleCheck
//      (
//        result( 0).Third() <= 0.8,
//        true
//      );
      return 0;
    }

    static const ExampleClass::EnumType ExampleChemistryConformationComparisonMolAlignByParts_Instance;

  }; //end ExampleChemistryMolAlignByParts

  const ExampleClass::EnumType ExampleChemistryMolAlignByParts::ExampleChemistryConformationComparisonMolAlignByParts_Instance
  (
    GetExamples().AddEnum( ExampleChemistryMolAlignByParts())
  );

} // namespace bcl
