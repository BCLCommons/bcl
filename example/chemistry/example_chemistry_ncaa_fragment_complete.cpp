// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
// (c) This file is part of the BCL software suite and is made available under license.
// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
// (c) For academic and non-profit users: 
// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
// (c) For commercial users: 
// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
// (c) For BCL developers at Vanderbilt University: 
// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
// (c)
// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the 
// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
// (c) Vanderbilt University. 
// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org 
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)

// include example header

#include "example.h"

// include the header of the class which this example is for
#include "chemistry/bcl_chemistry_ncaa_fragment_complete.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_factory.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "sdf/bcl_sdf.h"

using bcl::sdf::FragmentFactory;

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_ncaa_fragment_complete.cpp
  //! @details Tests ChemistryNCAAFragmentComplete class which contains small molecule configuration data
  //!
  //! @author vuot2, brownbp1, mendenjl
  //! @date
  //! @remarks status complete
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryNCAAFragmentComplete :
    public ExampleInterface
  {
  public:

    ExampleChemistryNCAAFragmentComplete *Clone() const
    {
      return new ExampleChemistryNCAAFragmentComplete( *this);
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

      // construct from conformation interface
      io::IFStream input_sdf;
      const std::string hexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, hexane_filename);

      // load information into small_mol_conformation
      chemistry::NCAAFragmentComplete small_mol_conformation( FragmentFactory::MakeFragment(input_sdf, sdf::e_Maintain));
      // close the input stream
      io::File::CloseClearFStream( input_sdf);

    /////////////////
    // data access //
    /////////////////

      // check GetAtomsIterator
      size_t atom_index = 0;
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface>
        itr( small_mol_conformation.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        if( itr->GetAtomType() != chemistry::GetAtomTypes().C_TeTeTeTe)
        {
          break;
        }
        ++atom_index;
      }
      // there should be 18 atoms of C_TeTeTeTe type
      BCL_ExampleIndirectCheck( atom_index, 6, "Conformational iterator succeeds to iterate over underlying atom vector");

      // test SetName and GetName
      small_mol_conformation.SetName( "hexane");
      BCL_ExampleCheck( small_mol_conformation.GetName(), "hexane");

      // test GetBonds
      BCL_ExampleCheck( small_mol_conformation.GetNumberBonds(), 18);

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

  const ExampleClass::EnumType ExampleChemistryNCAAFragmentComplete::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryNCAAFragmentComplete())
  );

} // namespace bcl
