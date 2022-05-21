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
#include "descriptor/bcl_descriptor_atom_lone_pair_en.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_atom_lone_pair_en.cpp
  //!
  //! @author kothiwsk, mendenjl, geanesar
  //! @date Dec 14, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAtomLonePairEN :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAtomLonePairEN *Clone() const
    {
      return new ExampleDescriptorAtomLonePairEN( *this);
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

      // default constructor
      descriptor::AtomLonePairEN lp_en;

      // copy constructor
      descriptor::AtomLonePairEN lp_en_copy( lp_en);

      // make an atom property that also gets lp_en
      descriptor::CheminfoProperty lp_en_from_atom_properties( descriptor::GetCheminfoProperties().calc_LonePairEN);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( lp_en.GetAlias(), "Atom_LonePairEN");

      BCL_ExampleCheck( lp_en.GetString(), "Atom_LonePairEN");

    ///////////////
    // operators //
    ///////////////

      std::string filename( "CSD_first10mols_inp.sdf");

      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, filename));
      // read in ensemble
      chemistry::FragmentEnsemble ensemble( input);
      // close stream
      io::File::CloseClearFStream( input);

      for
      (
        storage::List< chemistry::FragmentComplete>::iterator
          itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        linal::Vector< float> lone_pair_electronegativities
        (
          lp_en_from_atom_properties->CollectValuesOnEachElementOfObject( *itr)
        );

        BCL_ExampleAssert( lone_pair_electronegativities.GetSize(), itr->GetNumberAtoms());
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( lp_en, lp_en_copy),
        true,
        "SigmaCharge I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAtomLonePairEN

  const ExampleClass::EnumType ExampleDescriptorAtomLonePairEN::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAtomLonePairEN())
  );

} // namespace bcl
