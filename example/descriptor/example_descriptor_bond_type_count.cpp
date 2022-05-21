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
#include "descriptor/bcl_descriptor_bond_type_count.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_bond_type_count.cpp
  //!
  //! @author mendenjl
  //! @date Apr 10, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorBondTypeCount :
    public ExampleInterface
  {
  public:

    ExampleDescriptorBondTypeCount *Clone() const
    {
      return new ExampleDescriptorBondTypeCount( *this);
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

      // property to count aromatic bonds
      descriptor::BondTypeCount nbondtype_isaromatic( chemistry::ConfigurationalBondTypeData::e_IsAromatic, 1);

      // property to count all bonds
      descriptor::BondTypeCount nbondtype_identity( chemistry::ConfigurationalBondTypeData::e_Identity, 1);

      // copy constructor
      descriptor::BondTypeCount nbondtype_copy( nbondtype_isaromatic);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( nbondtype_identity.GetAlias(), "BondTypeCount");
      BCL_ExampleCheck( nbondtype_isaromatic.GetAlias(), "BondTypeCount");

      BCL_ExampleCheck( nbondtype_identity.GetString(), "BondTypeCount(property=Identity,value=1)");

    ///////////////
    // operators //
    ///////////////

      // map from filename to # of expected bonds of a bondtype
      storage::Map< std::string, size_t> filenames_to_expected_n_aromatic_bonds;

      // diazepam has 12 aromatic bonds
      filenames_to_expected_n_aromatic_bonds[ "diazepam.sdf"] = 24;
      // taxol has 18 aromatic bonds
      filenames_to_expected_n_aromatic_bonds[ "taxol.sdf"] = 36;
      // hexane has no aromatic bonds
      filenames_to_expected_n_aromatic_bonds[ "hexane.sdf"] = 0;

      for
      (
        storage::Map< std::string, size_t>::const_iterator
          itr_files( filenames_to_expected_n_aromatic_bonds.Begin()),
          itr_files_end( filenames_to_expected_n_aromatic_bonds.End());
        itr_files != itr_files_end;
        ++itr_files
      )
      {
        // create input stream for reading a smallmolecule ensemble
        io::IFStream input;
        BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, itr_files->first));
        // read in ensemble
        chemistry::FragmentEnsemble ensemble( input);
        // close stream
        io::File::CloseClearFStream( input);

        const chemistry::FragmentComplete &mol( ensemble.GetMolecules().FirstElement());

        BCL_ExampleIndirectCheck
        (
          nbondtype_isaromatic.SumOverObject( mol).First(),
          float( itr_files->second),
          "Counting aromatic bonds in " + itr_files->first
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( nbondtype_identity, nbondtype_isaromatic),
        true,
        "NBondType I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorBondTypeCount

  const ExampleClass::EnumType ExampleDescriptorBondTypeCount::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorBondTypeCount())
  );

} // namespace bcl
