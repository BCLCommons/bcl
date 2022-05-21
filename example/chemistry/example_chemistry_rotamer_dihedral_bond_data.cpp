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
#include "chemistry/bcl_chemistry_rotamer_dihedral_bond_data.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_rotamer_dihedral_bond_data.cpp
  //! @details Tests ChemistryRotamerDihedralBondData class which contains small molecule configuration data
  //!
  //! @author kothiwsk
  //! @date
  //! @remarks status complete
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryRotamerDihedralBondData :
    public ExampleInterface
  {
  public:

    ExampleChemistryRotamerDihedralBondData *Clone() const
    {
      return new ExampleChemistryRotamerDihedralBondData( *this);
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

      // create input stream for reading molecule of interest which needs to generated
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));
      // read in molecule
      chemistry::FragmentEnsemble molecule( input, sdf::e_Remove);
      // close stream
      io::File::CloseClearFStream( input);

      // create a tree searcher
      chemistry::SearchFragmentLibraryFromTree search_lib
      (
        ( chemistry::RotamerLibraryFile( chemistry::RotamerLibraryInterface::GetDefault()))
      );

      // get molecule priority
      chemistry::PriorityDihedralAngles molecule_priority;
      storage::Pair< storage::Vector< double>, storage::Vector< storage::VectorND< 4, size_t> > > priority_information
      (
        molecule_priority( molecule.GetMolecules().FirstElement())
      );
      // get isomorphism between of fragments with molecule , for fragments contained in molecule
      util::ShPtrVector< chemistry::SmallMoleculeFragmentIsomorphism> fragment_isomorphisms
      (
        search_lib.FindFragmentsOfMolecule( molecule.GetMolecules().FirstElement())
      );

      // get the rotamer dihedral bond data for the molecule
      const util::ShPtrVector< chemistry::RotamerDihedralBondData> bond_mapping
      (
        chemistry::SmallMoleculeFragmentMapping().MapFragmentIsomorphisms
        (
          molecule.GetMolecules().FirstElement(), fragment_isomorphisms
        )
      );
      io::OFStream os;
      const std::string fname( AddExampleOutputPathToFilename( molecule_priority, "bond_dihedral_data.txt"));
      BCL_ExampleMustOpenOutputFile( os, fname);
      io::Serialize::Write( bond_mapping, os);
      io::File::CloseClearFStream( os);

      if
      (
        BCL_ExampleCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance( fname, fname + ".correct.bz2", 0.01),
          true
        )
      )
      {
        remove( fname.c_str());
      }

    /////////////////
    // data access //
    /////////////////

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

  }; //end ExampleSmallMolecule

  const ExampleClass::EnumType ExampleChemistryRotamerDihedralBondData::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryRotamerDihedralBondData())
  );

} // namespace bcl
